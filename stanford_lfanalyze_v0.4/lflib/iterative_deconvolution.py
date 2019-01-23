# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

from lflib.lightfield import LightField
from lflib.imageio import save_image

import numpy as np
import time

from scipy.sparse.linalg.interface import LinearOperator
from scipy.sparse.dia import dia_matrix

from pylflib import compute_long_object_compensation

# -------- Variance stabilizing transforms -----------                                                                                                                       

def anscombe(x):
   return 2.0*np.sqrt(x + 3.0/8.0)

def inverse_anscombe(z):
   return (z/2.0)**2 - 3.0/8.0

def generalized_anscombe(x,mu,sigma,gain=1.0):
   return (2.0/gain)*np.sqrt(gain*x + (gain**2)*3.0/8.0 + sigma**2 - gain*mu)

def inverse_generalized_anscombe(z,mu,sigma,gain=1.0):
   return (gain*z/2.0)**2 - (gain**2)*3.0/8.0 - sigma**2 + gain * mu

# ------------------------------------------------------

# TODO:
#
# Not sure if this is necessary, but it's a handy snippet of code.
#
#    # Scale A by dividing by its largest eigenvalue.  This keeps the algorithm stable!
#    print '\t--> Estimating largest eigenvalue of A using power iteration.'
#    (A,b,scale_factor) = ScaleA(A,b)
#    print '\t    Scaled A by the its largest eigenvalue: ',  1.0/np.square(scale_factor)
#
def ScaleA(A,b):
    '''
    Returns: [A,b,scale_factor]

    Scales mu, A and f so that the largest eigenvalue of A.T*A is 1 and the
    new problem
 
    min sum_i (\|w_i|\ + beta/2 \|D_i u - w_i\|^2) + mu/2 \|Au - b\|^2
    
    is equivalent to the old one.  
    '''
    tol = .05;
    max_iterations = 10;

    # PORT : Handle real complex matrices here!
    #
    # if not isreal(A(rand(n,1),1)):
    #     eopts.isreal = false;

    # Compute largest eigenvalue by power iteration
    x = np.random.rand(A.shape[1])
    iteration = 0
    norm_z = tol + 1;
    while (iteration < max_iterations) and (norm_z > tol):
        z = A.rmatvec(A.matvec(x))
        norm_z = np.linalg.norm(z)
        print norm_z
        x = z / norm_z
        iteration += 1

    if np.real(norm_z) > 1 + 1e-10:
        b = b/np.sqrt(norm_z);
        A_operator = ScaledLinearOperator(A, 1.0/np.sqrt(norm_z))
        A = LinearOperator(A.shape,	  	
                           matvec=A_operator.matvec,
                           rmatvec=A_operator.rmatvec,
                           dtype='float')

    return (A, b, 1.0/np.sqrt(norm_z))

class ScaledLinearOperator(object):
    def __init__(self, A, scale_factor):
        self.A = A
        self.scale_factor = scale_factor

    def matvec(self, x ):
        return self.A.matvec(x) * self.scale_factor

    def rmatvec(self, x ):
        return self.A.rmatvec(x) * self.scale_factor

class LightFieldOperator(object):
    def __init__(self, sirt, db):
        self.sirt = sirt
        self.db = db
        self.left_preconditioner = None

    def matvec(self, vol_vec ):
        vol = np.reshape(vol_vec, (self.db.ny, self.db.nx, self.db.nz))
        im = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        b = np.reshape(im, (im.shape[0]*im.shape[1]))
        if self.left_preconditioner == None:
            return b
        else:
            return self.left_preconditioner * b

    def rmatvec(self, lf_vec ):
        if self.left_preconditioner == None:
            lf = np.reshape(lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        else:
            lf = np.reshape(self.left_preconditioner*lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                          representation = LightField.TILED_SUBAPERTURE))
        return np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))

    def as_linear_operator(self, nrays, nvoxels):
        return LinearOperator((nrays, nvoxels),
                              matvec=self.matvec,
                              rmatvec=self.rmatvec,
                              dtype='float')


class NormalEquationLightFieldOperator(object):
    def __init__(self, sirt, db):
        self.sirt = sirt
        self.db = db

    def matvec(self, vol_vec ):
        vol = np.reshape(vol_vec.astype(np.float32), (self.db.ny, self.db.nx, self.db.nz))
        lf = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                          representation = LightField.TILED_SUBAPERTURE))
        return np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))

class RegularizedNormalEquationLightFieldOperator(object):
    def __init__(self, sirt, db, regularization_lambda, left_preconditioner = None):
        self.sirt = sirt
        self.db = db
        self.regularization_lambda = regularization_lambda
        self.left_preconditioner = None

    def matvec(self, vol_vec ):
        input_vol = np.reshape(vol_vec.astype(np.float32), (self.db.ny, self.db.nx, self.db.nz))
        lf = self.sirt.project(input_vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        if self.left_preconditioner != None:
           lf *= np.square(left_preconditioner)
        output_vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv,
                                                      self.db.ns, self.db.nt,
                                                      representation = LightField.TILED_SUBAPERTURE))

        # L2-Norm Regularization
        output_vol += self.regularization_lambda * self.regularization_lambda * input_vol

        # Laplacian Regularization
        # lap_lambda = 15  # FIXME: HARD CODED FOR NOW!
        # nx = self.db.nx
        # ny = self.db.ny
        # nz = self.db.nz
        # lapvol = np.copy(input_vol)
        # lapvol[0:ny-1,:,:] -= 1/6.0 * input_vol[1:ny  ,:,:]
        # lapvol[1:ny  ,:,:] -= 1/6.0 * input_vol[0:ny-1,:,:]
        # lapvol[:,0:nx-1,:] -= 1/6.0 * input_vol[:,1:nx  ,:]
        # lapvol[:,1:nx  ,:] -= 1/6.0 * input_vol[:,0:nx-1,:]
        # lapvol[:,:,0:nz-1] -= 1/6.0 * input_vol[:,:,1:nz  ]
        # lapvol[:,:,1:nz  ] -= 1/6.0 * input_vol[:,:,0:nz-1]
        # # Zero out laplacian around the edges.
        # lapvol[0,:,:] = 0.0;
        # lapvol[:,0,:] = 0.0;
        # lapvol[:,:,0] = 0.0;
        # lapvol[ny-1,:,:] = 0.0;
        # lapvol[:,nx-1,:] = 0.0;
        # lapvol[:,:,nz-1] = 0.0;
        # output_vol += lap_lambda * lap_lambda * lapvol
        
        return np.reshape(output_vol, (self.db.nx * self.db.ny * self.db.nz))

class CgIterationFunctor:
    def __init__(self, algorithm_name, linear_operator, b, nvoxels, nrays):
        self.algorithm_name = algorithm_name
        self.iterations = 0
        self.prev_x = None
        self.last_time = 0
        self.linear_operator = linear_operator
        self.b = b
        self.nvoxels = nvoxels
        self.nrays = nrays

    def iter_callback(self, x):
        toc = time.time()
        
        if self.prev_x != None:
            print '\t--> [ CG Iteration %d   (%0.2f seconds) ] ' % (self.iterations,
                                                                   toc - self.last_time)
            if self.linear_operator != None:
                residual_norm = np.linalg.norm(self.b - self.linear_operator.matvec(x)) / self.nrays
                print '\t      Residual Norm: %0.4g' % (residual_norm)
                update_norm = np.linalg.norm(self.prev_x-x) / self.nvoxels
                print '\t        Update Norm: %0.4g' % (update_norm)
        
        self.last_time = toc
        self.prev_x = np.copy(x)
        self.iterations += 1

# ----------------------------------------------------------------------------------------
#                         EXPERIMENTAL DECONVOLUTION METHODS
# ----------------------------------------------------------------------------------------

def conjugate_gradient_reconstruction(lfcal, lightfield, alpha,
                                      convergence_threshold, max_iterations, regularization_lambda, 
                                      disable_gpu = False, gpu_id = 0):

    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db

    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    lfproj.set_premultiplier(lfcal.radiometric_correction)

    # Form the b matrix 
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    # Apply the preconditioner
    preconditioner = None
    #    preconditioner = 1.0 / np.sqrt(7.0 + b)  # Read noise variance of 7 for now...
    if preconditioner != None:
       b *= preconditioner

    # Uncomment to swap between zero initial volume and focal stack
    # volume.  (Zero volume start seems to converge faster, though).
    vol_vec = np.zeros((db.ny*db.nx*db.nz), dtype=np.float32)
    #vol_vec = np.reshape(lfproj.backproject(lightfield), (db.nx * db.ny * db.nz))

    # Conjugate gradient requires a square A matrix, so we solve the
    # normal equations below, rather than the original problem Ax = b.
    A_operator = LightFieldOperator(lfproj, db)
    AtA_operator = RegularizedNormalEquationLightFieldOperator(lfproj, db, regularization_lambda, preconditioner)
    A = LinearOperator((db.nt*db.nv*db.ns*db.nu, db.nx*db.ny*db.nz),
                       matvec=A_operator.matvec, rmatvec=A_operator, dtype='float')
    AtA = LinearOperator((db.nx*db.ny*db.nz, db.nx*db.ny*db.nz),
                         matvec=AtA_operator.matvec, rmatvec=AtA_operator.matvec, dtype='float')
    At_b = A_operator.rmatvec(b)

    print 'Calling Conjugate Gradient solver...'
    from scipy.sparse.linalg import cg

    nrays = db.ns*db.nt*db.nu*db.nv
    nvoxels = db.nx*db.ny*db.nz
    iter_func = CgIterationFunctor("CG", A, b, nvoxels, nrays)
    (vol_vec, info) = cg(AtA, At_b, x0 = vol_vec, tol = convergence_threshold, maxiter = max_iterations, callback=iter_func.iter_callback)

    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    vol = np.reshape(vol_vec, (db.ny, db.nx, db.nz))
    min_val = vol[db.supersample_factor:-db.supersample_factor,
                  db.supersample_factor:-db.supersample_factor, :].min()
    print '\t--> Replacing border values with min value: ', min_val
    vol[0:db.supersample_factor, :, :] = min_val
    vol[-db.supersample_factor:, :, :] = min_val
    vol[:, 0:db.supersample_factor, :] = min_val
    vol[:, -db.supersample_factor:, :] = min_val
    return vol.astype(np.float32)

# ------------------------------- ADMM TOTAL VARIATION SOLVER ----------------------------------

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
    
class AugmentedLightFieldOperator(object):
    def __init__(self, sirt, db, rho, structure_matrix):
        self.sirt = sirt
        self.db = db
        self.rho = rho
        self.structure_matrix = structure_matrix

    def matvec(self, vol_vec ):

        # Compute A*x
        vol = np.reshape(vol_vec, (self.db.ny, self.db.nx, self.db.nz))
        im = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        im_vec = np.reshape(im, (im.shape[0]*im.shape[1]))

        # Add the L2-Norm regularization term
        if self.structure_matrix != None:
            reg_vec = np.sqrt(self.rho) * self.structure_matrix * vol_vec
        else:
            reg_vec = np.sqrt(self.rho) * vol_vec

        return np.concatenate((im_vec, reg_vec), axis=0)

    def rmatvec(self, vec ):

        # Compute transpose(A)*x
        lf_vec_len = (self.db.ns*self.db.nt*self.db.nu*self.db.nv)
        lf_vec = vec[0:lf_vec_len]
        lf = np.reshape(lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                               representation = LightField.TILED_SUBAPERTURE))
        vol_vec = np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))

        # Compute rho * reg_vec
        if self.structure_matrix != None:
            reg_vec = np.sqrt(self.rho) * self.structure_matrix.T * vec[lf_vec_len:]
        else:
            reg_vec = np.sqrt(self.rho) * vec[lf_vec_len:]

        return vol_vec + reg_vec


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
    #    lfproj.set_premultiplier(lfcal.radiometric_correction)

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

    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    min_val = vol[db.supersample_factor:-db.supersample_factor,
                  db.supersample_factor:-db.supersample_factor, :].min()
    print '\t--> Replacing border values with min value: ', min_val
    vol[0:db.supersample_factor, :, :] = min_val
    vol[-db.supersample_factor:, :, :] = min_val
    vol[:, 0:db.supersample_factor, :] = min_val
    vol[:, -db.supersample_factor:, :] = min_val

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

# ------------------------------- TVAL3 SOLVER ---------------------------------
def tval3_reconstruction(db, lightfield, alpha,
                         convergence_threshold, max_iterations,
                         lambda_tv, 
                         initial_volume = None, debug = False,
                         disable_gpu = False, gpu_id = 0):

    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    nu = db.nu
    nv = db.nv
    ns = db.ns
    nt = db.nt

    # Produce the right hand side
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    # Create the linear operator
    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = LightFieldOperator(lfproj, db)
    A = LinearOperator((nrays, nvoxels),	  	
                       matvec=A_operator.matvec,
                       rmatvec=A_operator.rmatvec,
                       dtype='float')

    print 'Calling TVAL3 total variation solver...'

    ## Run TVAL3
    from PyTVAL3 import TVAL3
    opts = {}
    opts['mu'] = np.power(2,8);
    opts['beta'] = np.power(2,5);
    opts['tol'] = 1E-3;
    opts['maxit'] = 300;
    opts['TVnorm'] = 1;
    opts['nonneg'] = True;
    opts['TVL2'] = True;
    opts['disp'] = False;

    t = time.time();
    (U, out) = TVAL3(A,b,(db.ny,db.nx,db.nz),opts);
    print U.max()
    print U.min()
    save_image('test.tif', U, dtype=np.float32)
    t = time.time() - t;

        # toc = time.time()
        # print '\t--> [ ADMM total variation iteration %d,  rho = %0.2f   (%0.2f seconds) ] ' % (i, rho, toc-tic)
        # print '\t       Primal norm: %0.4g\t    (eps_pri  : %0.4g)' % (r_norm, eps_pri)
        # print '\t         Dual norm: %0.4g\t    (eps_dual : %0.4g)' % (s_norm, eps_dual)

        # # ---- DEBUG ----
        # A_operator_debug = LightFieldOperator(lfproj, db)
        # error = A_operator_debug.matvec(x) - b
        # nrays = db.ns*db.nt*db.nu*db.nv
        # residual_norm = np.linalg.norm(error) / nrays
        # print '\t     Residual Norm: %0.4g  ' % (residual_norm)
        # # ---------------

    return U.astype(np.float32)




# ------------------------------- LSQR SOLVER ----------------------------------

def lsqr_reconstruction(lfcal, lightfield, alpha,
                        convergence_threshold, max_iterations,
                        regularization_lambda,
                        initial_volume = None, debug = False,
                        disable_gpu = False, gpu_id = 0):

    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db

    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    lfproj.set_premultiplier(lfcal.radiometric_correction)


    # # DEBUG - MANUAL RADIOMETRY FOR TRYING VARIOUS STUFF OUT
    # lightfield_im = lightfield.asimage(representation = LightField.TILED_LENSLET)
    
    # vol_ones = db.ones_volume()
    # lf_ones = lfproj.project(vol_ones)
    # vol_weights = lfproj.backproject(lf_ones)
    # ideal_lf = lf_ones.asimage(representation = LightField.TILED_LENSLET)

    # rectified_radiometry = lfcal.rectified_radiometry

    # # # radiometric_correction = (rectified_radiometry) / (ideal_lf + 1e-16)
    # radiometric_correction = 1.0 / (rectified_radiometry)
    # radiometric_correction[np.nonzero(rectified_radiometry == 0)] = 0.0   # Prevent NaNs!
    # radiometric_correction[np.nonzero(ideal_lf == 0)] = 0.0   # Prevent NaNs!
    # lightfield_im *= radiometric_correction
    
    # lightfield = LightField(lightfield_im, db.nu, db.nv, db.ns, db.nt,
    #                         representation = LightField.TILED_LENSLET)
    # # # /DEBUG

    # Generate the right hand side vector
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = LightFieldOperator(lfproj, db)
    A = A_operator.as_linear_operator(nrays, nvoxels)

    print 'Calling LSQR solver...'
    from scipy.sparse.linalg import lsqr

    result = lsqr(A, b, damp = regularization_lambda, iter_lim = max_iterations, show=True)
    vol = np.reshape(result[0], (db.ny, db.nx, db.nz))

    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    min_val = vol[db.supersample_factor:-db.supersample_factor,
                  db.supersample_factor:-db.supersample_factor, :].min()
    print '\t--> Replacing border values with min value: ', min_val
    vol[0:db.supersample_factor, :, :] = min_val
    vol[-db.supersample_factor:, :, :] = min_val
    vol[:, 0:db.supersample_factor, :] = min_val
    vol[:, -db.supersample_factor:, :] = min_val

    return vol.astype(np.float32)


# ----------------------------------------------------------------------------------------
#                            SIRT ITERATIVE DECONVOLUTION
# ----------------------------------------------------------------------------------------



def sirt_reconstruction(lfcal, lightfield, alpha,
                        convergence_threshold, max_iterations,
                        regularization_lambda,
                        debug = False,
                        long_object = False,
                        disable_gpu = False, gpu_id = 0,save_errors=False,
                        debug_path = 'sirt_debug'):

    # Prefer wave optics model over geometric optics model
    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db
        
    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    lfproj.set_premultiplier(lfcal.radiometric_correction)
    
    # DEBUG - MANUAL RADIOMETRY FOR TRYING VARIOUS STUFF OUT

    #lightfield_im = lightfield.asimage(representation = LightField.TILED_LENSLET)
    #save_image("lambda.tif", lightfield_im, dtype=np.float32);
    
    #lightfield_im /= 10
    #lightfield = LightField(lightfield_im, db.nu, db.nv, db.ns, db.nt, representation = LightField.TILED_LENSLET)
    
    # vol_ones = db.ones_volume()
    # lf_ones = lfproj.project(vol_ones)
    # vol_weights = lfproj.backproject(lf_ones)
    # ideal_lf = lf_ones.asimage(representation = LightField.TILED_LENSLET)
    # # save_image("ideal_im.jpg", ideal_lf / ideal_lf.max()*255, dtype=np.uint8)
            
    # rectified_radiometry = lfcal.rectified_radiometry
    # self.radiometric_correction = self.rectified_radiometry / (self.ideal_lf + 1e-16)

    # /DEBUG




    # Create a linear operator for the optical model A.  This model
    # allows us to copmute A or A.T by calling its matvec() and
    # rmatvec() methods.
    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = LightFieldOperator(lfproj, db)
    A = A_operator.as_linear_operator(nrays, nvoxels)
    
    # Generate the b vector, which contains the observed lightfield;
    # and the initial volume x containing all zeros.
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, np.prod(im_subaperture.shape))
    x = np.zeros((nvoxels), dtype=np.float32)

    # Anscombe params
    ANSCOMBE = 0
    camera_gain = 2.0;
    readnoise_mean = 0.0;
    readnoise_sigma = 2.51;

    if ANSCOMBE == 1:
        print '\t--> Stabilizing variance using the anscombe transform'
        b = anscombe(b)
    elif ANSCOMBE == 2:
        print '\t--> Stabilizing variance using the generalized anscombe transform'
        b = generalized_anscombe(b, readnoise_mean, readnoise_sigma, camera_gain)
    else:
        print '\t--> Variance stabilization disabled.'



    # Model photon shot noise by setting the noise variance at every
    # sensor pixel to be equal to that pixels intenstity plus a
    # stationary read noise term.  In essence, this approximates the
    # photon shot noise + read noise as a non-stationary gaussian
    # distribution.  See this paper for more on this approximation:
    #
    # L. Mugnier et. al. MISTAL: a myopic edg-preserving image
    # restoration method, and application to astronomical
    # adaptive-optics-corrected long-exposure images.
    #

    # The Neo has 2e- RMS read noise, but I measured it to be 2.81 DN so we use that here.
    readnoise_variance = 7.84;  
    #photon_shotnoise_variance = b

    # DEBUG: Experimental code for low SNR zebrafish images
    #    
    # from lflib.imageio import load_image
    # variance_lf = LightField(load_image("lambda.tif"), db.nu, db.nv, db.ns, db.nt,
    #                          representation = LightField.TILED_LENSLET)
    # variance_photon = variance_lf.asimage(LightField.TILED_SUBAPERTURE)
    # photon_shotnoise_variance = np.reshape(variance_photon, np.prod(im_subaperture.shape)) 
    # A_operator.left_preconditioner = 1.0/np.square(photon_shotnoise_variance + readnoise_variance) 
    # A = A_operator.as_linear_operator(nrays, nvoxels)
    # b *= 1.0/np.square(photon_shotnoise_variance + readnoise_variance)

    # preconditioner = 1.0/np.square(b + readnoise_variance) 
    # A_operator.left_preconditioner = preconditioner
    # A = A_operator.as_linear_operator(nrays, nvoxels)
    # b *= preconditioner




    # Create the SIRT weight volume and light field.  These are
    # created by projecting a volume containing all ones (to create
    # the weight lightfield), and then back-projecting a light field
    # with all ones (to create the weight volume).
    #
    # x_weights = A.T * b_ones;    b_weights = A * x_ones
    import time

    tic = time.time()
    b_weights = A.matvec(np.ones((nvoxels), dtype=np.float32))
    print '\t--> Time for one forward projection: %0.2f seconds.' % (time.time()-tic)

    tic = time.time()
    x_weights = A.rmatvec(np.ones((nrays), dtype=np.float32))
    print '\t    Time for one back projection: %0.2f seconds.' % (time.time()-tic)

    # --------------------------------------------------------------------

    # Make sure that any zero valued weights are set to nonzero so
    # that they don't lead to division by zero below. We then
    # normalize the starting volume using the volume weights.
    min_bweight = b_weights[np.nonzero(b_weights != 0.0)].min()
    min_xweight = x_weights[np.nonzero(x_weights != 0.0)].min()
    b_weights[np.nonzero(b_weights < min_bweight)] = min_bweight;
    x_weights[np.nonzero(x_weights < min_xweight)] = min_xweight;    


    iteration_error = []
    for i in range(max_iterations):
        tic = time.time()

        # In each iteration, forward and backproject error from all views at once, then update the volume
        #
        # STEP 1: forward projection of volume to create sub-aperture images.
        b_hat = A.matvec(x)

        # DEBUGGING
        # if i == 1:
        #     b_debug = np.reshape(b_hat, (db.nt*db.nv, db.ns*db.nu))
        #     save_image("lf_" + str(i) + ".tif", b_debug);

        # STEP 2: Compute error between computed and observed sub-aperture images
        error = b - b_hat

        if i >= 1:
            error[np.nonzero(b_hat == 0.0)] = 0.0

        # Debug: save error images
        #
        #if i == max_iterations-1:
        error_test = np.reshape(inverse_anscombe(error), (db.nv*db.nt, db.nu*db.ns))
        save_image("error_" + str(i) + ".tif", error_test, dtype=np.float32);


        if i == 14:
            print 'DEBUG TIME'
            save_image("bhat_b_" + str(i) + ".tif", np.reshape(inverse_anscombe(b_hat),
                                                               (db.nt*db.nv, db.ns*db.nu)), dtype=np.float32);
            
        # collect the unweighted error in light field space
        if save_errors:
            iteration_error.append( error )

        # the reweighted error
        reweighted_error = error / b_weights

        # STEP 3: back-project error onto the volume
        if ANSCOMBE == 1:
            error_backprojected = A.rmatvec(inverse_anscombe(reweighted_error))
        elif ANSCOMBE == 2:
            error_backprojected = A.rmatvec(inverse_generalized_anscombe(reweighted_error,
                                                                         readnoise_mean,
                                                                         readnoise_sigma,
                                                                         camera_gain))
        else:
            error_backprojected = A.rmatvec(reweighted_error)

        # Graph Laplacian Regularization
        #
        # WARNING: This code has not been carefully checked and may not work!  -broxton
        # if regularization_lambda > 0.0:
        #     vol = np.reshape(x, (db.ny, db.nx, db.nz))
        #     lapvol = np.copy(vol)
        #     lapvol[0:ny-1,:,:] -= 1/6.0 * vol[1:ny  ,:,:]
        #     lapvol[1:ny  ,:,:] -= 1/6.0 * vol[0:ny-1,:,:]
        #     lapvol[:,0:nx-1,:] -= 1/6.0 * vol[:,1:nx  ,:]
        #     lapvol[:,1:nx  ,:] -= 1/6.0 * vol[:,0:nx-1,:]
        #     lapvol[:,:,0:nz-1] -= 1/6.0 * vol[:,:,1:nz  ]
        #     lapvol[:,:,1:nz  ] -= 1/6.0 * vol[:,:,0:nz-1]
        #     # Zero out laplacian around the edges.
        #     lapvol[0,:,:] = 0.0;
        #     lapvol[:,0,:] = 0.0;
        #     lapvol[:,:,0] = 0.0;
        #     lapvol[ny-1,:,:] = 0.0;
        #     lapvol[:,nx-1,:] = 0.0;
        #     lapvol[:,:,nz-1] = 0.0;
        #     lapvol_vec = np.reshape(lapvol, (db.ny*db.nx*db.nz))

        #     # Apply the reweighting and step size.  
        #     x_update = (error_backprojected / x_weights - regularization_lambda * lapvol_vec)
        # else:
        #     x_update = (error_backprojected / x_weights)

        # L2 Penalty (Tihkonov) Regularization
        #
        if regularization_lambda > 0.0:
            x_update = (error_backprojected / x_weights - regularization_lambda * x * x)
        else:
            x_update = error_backprojected / x_weights

        # Apply the update
        x += alpha * x_update

        # Debugging
        # if i == 1:
        #     vol_debug = np.reshape(x, (db.ny, db.nx, db.nz))
        #     save_image("vol_" + str(i) + ".tif", vol_debug)

        # Enforce non-negativity constraint
        x[x<0]=0

        # CHECK FOR CONVERGENCE
        #
        # normalize MSE using input LF
        nrays = db.ns*db.nt*db.nu*db.nv
        residual_norm = np.linalg.norm(error) / nrays

        # normalize MSE using input LF
        nvoxels = db.nx*db.ny*db.nz
        update_norm = np.linalg.norm(alpha * x_update) / nvoxels

        toc = time.time()
        print '\t--> [ SIRT Iteration %d   (%0.2f seconds) ] ' % (i, toc-tic)
        print '\t      Residual Norm: %0.4g' % (residual_norm)
        print '\t        Update Norm: %0.4g               (tol = %0.2e)  ' % (update_norm, convergence_threshold)

        # check if convergence criteria met
        if i > 0 and update_norm < convergence_threshold:
            break 

        # save out iteration errors
        if save_errors:
            np.savez( "iteration_errors", iteration_error )

    vol = np.reshape(x, (db.ny, db.nx, db.nz)) # Note that the default order for np.reshape in 'C' (row-major)

    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    min_val = vol[db.supersample_factor:-db.supersample_factor,
                  db.supersample_factor:-db.supersample_factor, :].min()
    print '\t--> Replacing border values with min value: ', min_val
    vol[0:db.supersample_factor, :, :] = min_val
    vol[-db.supersample_factor:, :, :] = min_val
    vol[:, 0:db.supersample_factor, :] = min_val
    vol[:, -db.supersample_factor:, :] = min_val
    return vol.astype(np.float32)

# ----------------------------------------------------------------------------------------
#                            AMP ITERATIVE DECONVOLUTION
# ----------------------------------------------------------------------------------------


def amp_reconstruction(lfcal, lightfield, alpha,
                       convergence_threshold, max_iterations,
                       regularization_lambda,
                       debug = False,
                       long_object = False,
                       disable_gpu = False, gpu_id = 0,save_errors=False,
                       debug_path = 'sart_debug',
                       wavelet_smoothing = True):

    # import wavelet functions if needed 
    if wavelet_smoothing:
        from lflib.wavelet3d import undecimated_3d, inverse_undecimated_3d, modify_wavelet_coefs, wavelet_threshold

    # Prefer wave optics model over geometric optics model
    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db
        
    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    lfproj.set_radiometry_correction(lfcal.radiometric_correction)

    nu = db.nu
    nv = db.nv
    ns = db.ns
    nt = db.nt

    # Generate the b vector, which contains the observed lightfield 
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    # Create a linear operator for the optical model A.  This model
    # allows us to copmute A or A.T by calling its matvec() and
    # rmatvec() methods.
    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = LightFieldOperator(lfproj, db)
    A = A_operator.as_linear_operator(nrays, nvoxels)

    # Model photon shot noise by setting the noise variance at every
    # sensor pixel to be equal to that pixels intenstity.  This should
    # be true if photon shot noise is the dominating noise term.
    #
    EPS = 1e-1        # Avoid dividing by zero!   This value works well on fish volumes, but maybe needs tuning?
    A_operator.left_preconditioner = 1.0/np.sqrt(b+EPS) 
    A = A_operator.as_linear_operator(nrays, nvoxels)
    b *= 1.0/np.sqrt(b+EPS)

    # --------------------------------------------------------------------
    
    if save_errors:
        iteration_error = []

    for i in range(max_iterations):
        tic = time.time()

        # In each iteration, forward and backproject error from all views at once, then update the volume
        #
        # STEP 1: forward projection of volume to create sub-aperture images. 
        if wavelet_smoothing:
            # A \Phi x
            b_hat = A.matvec(inverse_undecimated_3d(x))
        else:
            # Ax
            b_hat = A.matvec(x)

        # DEBUGGING
        # if i == 1:
        #     b_debug = np.reshape(b_hat, (db.nt*db.nv, db.ns*db.nu))
        #     save_image("lf_" + str(i) + ".tif", b_debug);

        # STEP 2: Compute error between computed and observed sub-aperture images
        error = b - b_hat

        # collect the unweighted error in light field space
        if save_errors:
            iteration_error.append( error )

        # add AMP adjustment
        adjusted_error = error + (1.0/delta)*np.sum(last_estimate)*last_error

        # STEP 3: back-project error onto the volume or wavelet coefficient space
        if wavelet_smoothing:
            error_backprojected = undecimated_3d( A.rmatvec(reweighted_error) )
        else:
            error_backprojected = A.rmatvec(reweighted_error)

        # Apply the update -- should we be using x for wavelet coefs and volume interchangably,
        # or should we use a for the wavelet coefs if wavelet_smoothing = True?
        if wavelet_smoothing:
           x = modify_wavelet_coefs(x, error_backprojected, scale = alpha) # modify wavelets in R
           x = wavelet_threshold(x) # threshold modified coefs
        else:
           x_update = (error_backprojected / x_weights)
           x += alpha * x_update
        # Debugging
        # if i == 1:
        #     vol_debug = np.reshape(x, (db.ny, db.nx, db.nz))
        #     save_image("vol_" + str(i) + ".tif", vol_debug)

        # Enforce non-negativity constraint
        if wavelet_smoothing:
            x = wavelet_pos_thresh(x)
        else:
            x[x<0]=0

        # CHECK FOR CONVERGENCE
        #
        # normalize MSE using input LF
        nrays = db.ns*db.nt*db.nu*db.nv
        residual_norm = np.linalg.norm(error) / nrays

        # normalize MSE using input LF
        nvoxels = db.nx*db.ny*db.nz
        update_norm = np.linalg.norm(alpha * x_update) / nvoxels

        toc = time.time()
        print '\t--> [ SIRT Iteration %d   (%0.2f seconds) ] ' % (i, toc-tic)
        print '\t      Residual Norm: %0.4g' % (residual_norm)
        print '\t        Update Norm: %0.4g               (tol = %0.2e)  ' % (update_norm, convergence_threshold)

        # check if convergence criteria met
        if i > 0 and update_norm < convergence_threshold:
            break 

        # save out iteration errors
        if save_errors:
            np.savez( "iteration_errors", iteration_error )

    vol = np.reshape(x, (db.ny, db.nx, db.nz))

    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    vol[0:db.supersample_factor, :, :] = 0.0
    vol[-db.supersample_factor:, :, :] = 0.0
    vol[:, 0:db.supersample_factor, :] = 0.0
    vol[:, -db.supersample_factor:, :] = 0.0
    
    return vol.astype(np.float32)


