# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2013 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
import time

from scipy.sparse.linalg.interface import LinearOperator

from lflib.lightfield import LightField
from lflib.imageio import save_image
from lflib.linear_operators import LightFieldOperator, AugmentedLightFieldOperator

# ----------------------------------------------------------------------------------------
#                               ADMM ITERATIVE SOLVERS
# ----------------------------------------------------------------------------------------

# ---------------------------- ADMM TOTAL VARIATION SOLVER -------------------------------

def l1_shrinkage(x, kappa):
    return np.maximum( 0, x - kappa ) - np.maximum( 0, -x - kappa );

def sparse_total_variation_matrix(nx, ny, nz):
    '''
    Builds a sparse, square matrix that computes the total variation
    of a vector x.  The total variation matrix is defined as:
    
            n   for i==j where n is the number of voxels adjacent to x_i
    F_ij = -1   for i \neq j but adjacent to j
            0   otherwise
    '''
    nvoxels = nx*ny*nz

    # DEBUGGING: ADD A TINY BIT OF LASSO REGULARIZATION TO STABILIZE THE RECONSTRUCTION
    LASSO_REG = 0.0  # 0.1 seems to work well on 1x volumes

    from scipy.sparse import coo_matrix

    y_coords = np.reshape(np.tile(np.arange(ny, dtype=np.int32), (nx*nz, 1)), (nvoxels), order='f')
    x_coords = np.reshape(np.tile(np.reshape(np.tile(np.arange(nx, dtype=np.int32), (nz, 1)),
                                             (nz*nx), order='f'), (ny, 1)), (nx*ny*nz))
    z_coords = np.tile(np.arange(nz, dtype=np.int32), (nx*ny))

    # Linear index into coordinates involves combining 3D coordinates in (y, x, z) order (for now).
    diag_coords = y_coords*nx*nz + x_coords*nz + z_coords

    # Form the z+1 difference entries
    valid_idxs = np.nonzero(np.logical_and(z_coords-1 >= 0,  z_coords+1 < nz))
    diff_coords = y_coords*nx*nz + x_coords*nz + (z_coords+1)
    Dz = coo_matrix((np.ones(len(valid_idxs[0]))*1.0, (diag_coords[valid_idxs], diff_coords[valid_idxs])),
                     shape = (nvoxels, nvoxels), dtype=np.float32)

    # Form the z-1 difference entries
    diff_coords = y_coords*nx*nz + x_coords*nz + (z_coords-1)
    Dz = Dz + coo_matrix((np.ones(len(valid_idxs[0]))*-1.0, (diag_coords[valid_idxs], diff_coords[valid_idxs])),
                         shape = (nvoxels, nvoxels), dtype=np.float32)

    # Form the x+1 difference entries
    valid_idxs = np.nonzero(np.logical_and(x_coords-1 >= 0, x_coords+1 < nx))
    diff_coords = y_coords*nx*nz + (x_coords+1)*nz + z_coords
    Dx = coo_matrix((np.ones(len(valid_idxs[0]))*1.0, (diag_coords[valid_idxs], diff_coords[valid_idxs])),
                    shape = (nvoxels, nvoxels), dtype=np.float32)

    # Form the x-1 difference entries
    diff_coords = y_coords*nx*nz + (x_coords-1)*nz + z_coords
    Dx = Dx + coo_matrix((np.ones(len(valid_idxs[0]))*-1.0, (diag_coords[valid_idxs], diff_coords[valid_idxs])),
                       shape = (nvoxels, nvoxels), dtype=np.float32)

    # Form the y+1 difference entries
    valid_idxs = np.nonzero(np.logical_and(y_coords-1 >= 0, y_coords+1 < ny))
    diff_coords = (y_coords+1)*nx*nz + x_coords*nz + z_coords
    Dy = coo_matrix((np.ones(len(valid_idxs[0]))*1.0, (diag_coords[valid_idxs], diff_coords[valid_idxs])),
                       shape = (nvoxels, nvoxels), dtype=np.float32)

    # Form the y-1 difference entries
    diff_coords = (y_coords-1)*nx*nz + x_coords*nz + z_coords
    Dy = Dy + coo_matrix((np.ones(len(valid_idxs[0]))*-1.0, (diag_coords[valid_idxs], diff_coords[valid_idxs])),
                       shape = (nvoxels, nvoxels), dtype=np.float32)
    # Tests
    print 'Testing F matrix...'
    try:
        ones = np.ones((F.shape[0],))
        Finner1 = (F*ones).reshape(ny,nx,nz)[3:-3,3:-3,3:-3]
        print "\t--> Sum of F inner product with ones, cropped:", np.sum(Finner1)
        assert(np.sum(Finner1) == 0.0)
    except:
        print "\t--> First differences matrix has non-zero inner product with vector of ones!!!"
    try:
        ones = np.ones((F.shape[0],))
        Finner1 = (F.T*ones).reshape(ny,nx,nz)[3:-3,3:-3,3:-3]
        print "\t--> Sum of F.T inner product with ones, cropped:", np.sum(Finner1)
        assert(np.sum(Finner1) == 0.0)
    except:
        print "\t--> First differences matrix transpose has non-zero inner product with vector of ones!!!"

    # TODO: Think about relative weight of Z edges vs. X/Y edges.
    import scipy.sparse
    vs = scipy.sparse.vstack((Dy, Dx, Dz))
    return vs    

def admm_total_variation_reconstruction(lfcal, lightfield, alpha,
                                        convergence_threshold, max_iterations,
                                        lambda_tv, lambda_lasso, 
                                        initial_volume = None, debug = False,
                                        disable_gpu = False, gpu_id = 0):

    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db

    ABSTOL   = convergence_threshold        #  1e-4 seems to work well
    RELTOL   = 1e-2;
    ALPHA    = 1.5
    LSQR_ITERATIONS = 2
    ADMM_ITERATIONS = max_iterations

    # Penalty param hyperparameters
    rho = 0.2
    mu = 10.0
    tau_incr = 2.0
    tau_decr = 2.0

    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    nu = db.nu
    nv = db.nv
    ns = db.ns
    nt = db.nt

    # Produce the right hand side
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    F = sparse_total_variation_matrix(db.nx, db.ny, db.nz)
    if lambda_tv <= 0.0:
        print 'Error: you must specify a non-zero regularization value (lambda) for total variation.'
        exit(1)
    # if lambda_lasso > 0.0:
    #     print '\t--> Stabilizing solution with LASSO structure matrix (i.e. the identity)...'
    #     from scipy.sparse import eye
    #     nvoxels = db.nx * db.ny * db.nz
    #     F = F + lambda_lasso/lambda_tv * eye(nvoxels, nvoxels, dtype=np.float32)

    # Create the linear operator
    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = AugmentedLightFieldOperator(lfproj, db, rho, structure_matrix = F)
    A = LinearOperator((nrays + 3*nvoxels, nvoxels),	  	
                       matvec=A_operator.matvec,
                       rmatvec=A_operator.rmatvec,
                       dtype='float')

    print 'Calling ADMM/LSMR total variation solver...'
    from lflib.lsmr import lsmr

    z = np.zeros(3*(nvoxels), dtype=np.float32)
    u = np.zeros(3*(nvoxels), dtype=np.float32)
    x = np.zeros((nvoxels), dtype=np.float32)
    for i in range(ADMM_ITERATIONS):
        tic = time.time()

        warm_start = True
        if warm_start:
            # Solve the X update with LSQR (with warm start)
            rhs = np.concatenate(( b, np.sqrt(rho)*(z-u) ), axis=0)
            c = rhs - A.matvec(x)
            result = lsmr(A, c, damp = lambda_lasso, maxiter = LSQR_ITERATIONS, show=False)
            dx = result[0]
            x = x + dx
        else:
            # Solve the X update with LSQR (with warm start)
            rhs = np.concatenate(( b, np.sqrt(rho)*(z-u) ), axis=0)
            result = lsmr(A, rhs, damp = lambda_lasso, maxiter = LSQR_ITERATIONS, show=False)
            x = result[0]

        # Enforce non-negativity constraint
        x[x<0]=0

        # Update z & u
        z_old = z
        x_hat = ALPHA * F*x + (1 - ALPHA) * z_old;
        z = l1_shrinkage( x_hat + u, lambda_tv/rho )

        u = u + x_hat - z

        # Compute primary and dual residuals
        r_norm = np.linalg.norm(F*x - z)
        s_norm = np.linalg.norm(-rho * F.T * ( z - z_old ))

        # Check for convergence
        eps_pri  = np.sqrt(nvoxels)*ABSTOL + RELTOL*max(np.linalg.norm(F*x), np.linalg.norm(-z));
        eps_dual = np.sqrt(nvoxels)*ABSTOL + RELTOL*np.linalg.norm(rho*F.T*u);

        toc = time.time()
        print '\t--> [ ADMM total variation iteration %d,  rho = %0.2f   (%0.2f seconds) ] ' % (i, rho, toc-tic)
        print '\t       Primal norm: %0.4g\t    (eps_pri  : %0.4g)' % (r_norm, eps_pri)
        print '\t         Dual norm: %0.4g\t    (eps_dual : %0.4g)' % (s_norm, eps_dual)

        # ---- DEBUG ----
        A_operator_debug = LightFieldOperator(lfproj, db)
        error = A_operator_debug.matvec(x) - b
        nrays = db.ns*db.nt*db.nu*db.nv
        residual_norm = np.linalg.norm(error) / nrays
        print '\t     Residual Norm: %0.4g  ' % (residual_norm)
        # ---------------

        if r_norm < eps_pri and s_norm < eps_dual:
            print 'ADMM converged.'
            break;

        # Update the penalty parameter
        if r_norm > mu * s_norm:
            rho *= tau_incr
            u /= tau_incr    # The residual must be updated as well!
        elif s_norm > mu * r_norm:
            rho /= tau_decr
            u *= tau_decr    # The residual must be updated as well!
        A_operator.lambda_tv = rho


    # DEBUG: Save out gradient info
    print "\t--> Saving gradient volume (for debugging)"
    grad_vec = np.power(z[0:nvoxels],2) + np.power(z[nvoxels:2*nvoxels],2) + np.power(z[2*nvoxels:],2)
    grad_vol = np.reshape(grad_vec, (db.ny, db.nx, db.nz))
    save_image('gradvol.tif', grad_vol, dtype=np.float32)

    vol = np.reshape(x, (db.ny, db.nx, db.nz))
    return vol.astype(np.float32)

# ------------------------------- ADMM HUBER SOLVER ----------------------------------

def huber_shrinkage(x, kappa):
    return np.maximum(0, 1.0 - kappa/np.abs(x)) * x
    # result = vec.copy()
    # result[np.nonzero(vec > threshold)] -= threshold
    # result[np.nonzero(np.abs(vec) < threshold)] = 0.0
    # result[np.nonzero(vec < -threshold)] += threshold
    # return result

def admm_huber_reconstruction(lfcal, lightfield, alpha,
                              convergence_threshold, max_iterations,
                              regularization_lambda,
                              initial_volume = None, debug = False,
                              disable_gpu = False, gpu_id = 0):

    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db

    ABSTOL   = convergence_threshold        #  1e-4 seems to work well
    RELTOL   = 1e-2;
    ALPHA    = 1.8
    LSQR_ITERATIONS = 2
    ADMM_ITERATIONS = max_iterations

    # Penalty param hyperparameters
    rho = 0.01
    mu = 10.0
    tau_incr = 2.0
    tau_decr = 2.0

    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    nu = db.nu
    nv = db.nv
    ns = db.ns
    nt = db.nt

    # Produce the right hand side
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = LightFieldOperator(lfproj, db)
    A = LinearOperator((nrays, nvoxels),
                       matvec=A_operator.matvec,
                       rmatvec=A_operator.rmatvec,
                       dtype='float')

    print 'Calling ADMM/LSQR huber/l1 solver...'
    from lflib.lsmr import lsmr

    z = np.zeros((nrays), dtype=np.float32)
    u = np.zeros((nrays), dtype=np.float32)
    x = np.zeros((nvoxels), dtype=np.float32)
    for i in range(ADMM_ITERATIONS):
        tic = time.time()

        # Update RHS
        v = b + z - u
        result = lsmr(A, v, damp = np.sqrt(2*regularization_lambda/rho),
                      maxiter = LSQR_ITERATIONS, show=False)
        x = result[0]

        # Update z with relaxation
        zold = z
        Ax = A_operator.matvec(x)
        Ax_hat = alpha * Ax + (1-alpha) * (zold + b)
        tmp = Ax_hat - b + u

        # Huber loss
        z = rho/(rho+1.0) * tmp  + 1.0/(rho+1.0) * huber_shrinkage( tmp, 1.0+1.0/rho )


        # Enforce non-negativity constraint
        # z[z<0]=0


        u = u + Ax_hat - z - b

        # Compute primary and dual residuals
        r_norm = np.linalg.norm( Ax - z - b )
        s_norm = np.linalg.norm( -rho * A_operator.rmatvec( z - zold ) )

        # Check for convergence
        eps_pri  = np.sqrt(nvoxels)*ABSTOL + RELTOL*max(np.linalg.norm(Ax),
                                                        np.linalg.norm(-z),
                                                        np.linalg.norm(b));
        eps_dual = np.sqrt(nvoxels)*ABSTOL + RELTOL*np.linalg.norm(rho*u);

        toc = time.time()
        print '\t--> [ ADMM huber/l1 iteration %d,  rho = %0.2f   (%0.2f seconds) ] ' % (i, rho, toc-tic)
        print '\t       Primal norm: %0.4g\t    (eps_pri  : %0.4g)' % (r_norm, eps_pri)
        print '\t         Dual norm: %0.4g\t    (eps_dual : %0.4g)' % (s_norm, eps_dual)

        # ---- DEBUG ----
        A_operator_debug = LightFieldOperator(lfproj, db)
        error = A_operator_debug.matvec(x) - b
        nrays = db.ns*db.nt*db.nu*db.nv
        residual_norm = np.linalg.norm(error) / nrays
        print '\t     Residual Norm: %0.4g  ' % (residual_norm)
        # ---------------

        if r_norm < eps_pri and s_norm < eps_dual:
            print 'ADMM converged.'
            break;

        # Update the penalty parameter
        if r_norm > mu * s_norm:
            rho *= tau_incr
            u /= tau_incr    # The residual must be updated as well!
        elif s_norm > mu * r_norm:
            rho /= tau_decr
            u *= tau_decr    # The residual must be updated as well!

    vol = np.reshape(x, (db.ny, db.nx, db.nz))
    return vol.astype(np.float32)

#------------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#EOF
