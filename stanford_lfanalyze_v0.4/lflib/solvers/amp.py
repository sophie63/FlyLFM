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
from lflib.linear_operators import LightFieldOperator

# ----------------------------------------------------------------------------------------
#                APPROXIMATE MESSAGE PASSING (AMP) ITERATIVE SOLVER
# ----------------------------------------------------------------------------------------

def amp_reconstruction(lfcal, lightfield, 
                       alpha,
                       convergence_threshold, 
                       max_iterations,
                       sparsity,
                       delta,
                       debug = False,
                       long_object = False,
                       disable_gpu = False, 
                       gpu_id = 0,
                       save_errors=False,
                       debug_path = 'amp_debug',
                       multiscale_smoothing = False,
                       wavelet_type='la8',
                       stabilization=None,
                       standardize=True,
                       transform_type="undecimated_wavelet"):


    alpha = 1.0

    # import wavelet functions if needed 
    if multiscale_smoothing:
        from lflib.multiscale3d import ( multiscale_transform_3D, inverse_multiscale_transform_3D, 
                                         multiscale_coefficient_update, multiscale_threshold, 
                                         multiscale_coefficient_mean, generalized_anscombe,
                                         output_multiscale_coefs, anscombe, inverse_anscombe )

    # Prefer wave optics model over geometric optics model
    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db
        
    # get lightfield projection operator
    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, 
                                  gpu_id = gpu_id)
    lfproj.set_premultiplier(lfcal.radiometric_correction)

    # get lightfield dimensions
    nu = db.nu
    nv = db.nv
    ns = db.ns
    nt = db.nt

    # Generate the b vector, which contains the observed lightfield 
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    # Create a linear operator for the optical model A.  This model
    # allows us to compute A or A.T by calling its matvec() and
    # rmatvec() methods.
    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_operator = LightFieldOperator(lfproj, db)
    A = A_operator.as_linear_operator(nrays, nvoxels)

    if stabilization is not 'anscombe':
        # Model photon shot noise by setting the noise variance at every
        # sensor pixel to be equal to that pixels intensity.  This should
        # be true if photon shot noise is the dominating noise term.
        EPS = 1e-1        # Avoid dividing by zero!   
        # This value works well on fish volumes, but maybe needs tuning?

        A_operator.diagonalizer = 1.0/np.sqrt(b+EPS) 
        A = A_operator.as_linear_operator(nrays, nvoxels)
        b *= 1.0/np.sqrt(b+EPS)

    else: 
        print "\t--> Stabilizing variance with Anscombe transform."
        gain = 1.0
        mu = 0.0
        sigma = 10.0
        from lflib.multiscale3d import generalized_anscombe, anscombe, inverse_anscombe
#        b = generalized_anscombe(b, mu, sigma, gain)
        b = anscombe(b)
        print "Range b:", np.min(b), np.max(b)

    # Renormalize b to take values in [0,1] to improve numerical stability
    if standardize:
        lf_max = np.max(b)
        b /= lf_max
        print "Range b standardized:", np.min(b), np.max(b)

    # Initialize volume (or wavelet coefs)
    vol_shape = np.array([db.ny, db.nx, db.nz])
    if multiscale_smoothing:
        # construct set of zero-valued wavelet coefficients
        zero_vol = np.zeros((nvoxels), dtype=np.float32)
        J = int(np.ceil(np.log2(np.min(vol_shape))))
        print "J=",J
        x, vol_inds = multiscale_transform_3D( zero_vol, vol_shape=vol_shape, 
                                               wavelet_type=wavelet_type,
                                               transform_type=transform_type )

    else:
        # construct flattened zero volume 
        x = np.zeros((nvoxels), dtype=np.float32)

    # Initialize previous estimate average and residual vector
    last_estimate_avg = 0.0
    error = np.zeros((nrays), dtype=np.float32)

    # Create row and column sum vectors used to stabilize iterations.
    # These are created by projecting a volume containing all ones (to create
    # the weight lightfield), and then back-projecting a light field
    # with all ones (to create the weight volume).
    # They are used to reweight the estimated volume and errors in order 
    # to account for the fact that A and A.T have operator norms >1.
    #
    # N.B. -- there are probably better ways to precondition than this.
    x_ones = np.ones((nvoxels), dtype=np.float32)
    b_ones = np.ones((nrays), dtype=np.float32)

    # Compute weights:  x_weights = A.T * b_ones;    b_weights = A * x_ones
    import time
    print 'Computing forward weights...'
    tic = time.time()
    x_weights = A.rmatvec(b_ones)
    print '%0.2f seconds elapsed.' % (time.time()-tic)

    print 'Computing backward weights...'
    tic = time.time()
    b_weights = A.matvec(x_ones)
    print '%0.2f seconds elapsed.' % (time.time()-tic)
    
    # Make sure that any zero valued weights are set to nonzero so
    # that they don't lead to division by zero below. We then
    # normalize the starting volume using the volume weights.
    min_bweight = b_weights[np.nonzero(b_weights != 0.0)].min()
    min_xweight = x_weights[np.nonzero(x_weights != 0.0)].min()
    b_weights[np.nonzero(b_weights < min_bweight)] = min_bweight;
    x_weights[np.nonzero(x_weights < min_xweight)] = min_xweight;    

    #--------------------------------------------------------------------
    
    if save_errors:
        iteration_error = []

    for i in range(max_iterations):
        iteration_tic = time.time()

        # In each iteration, forward and backproject error 
        # from all views at once, then update the volume
        #
        # STEP 1: forward projection of volume to create sub-aperture images. 
        if multiscale_smoothing:
            # A \Phi x
            tic = time.time()
            if transform_type == "pyramidal_median":
                b_hat=A.matvec(np.asarray(inverse_multiscale_transform_3D(x,vol_inds=vol_inds,
                                transform_type=transform_type)).flatten() ) 
            else:
                inv_3d = np.asarray(inverse_multiscale_transform_3D(x,vol_inds=vol_inds,
                                                                    transform_type=transform_type))
                inv_3d = inv_3d.transpose((1,0,2)).flatten() # reshape for lflib 
                print "Range inv_3d:", np.min(inv_3d), np.max(inv_3d)
#                inv_3d[inv_3d < 0.0]=0.0 # impose nonnegativity
#                print "Max inv_3d:", np.max(inv_3d)
                b_hat = A.matvec( inv_3d ) 
                              
            b_hat /= b_weights
            print "\t--> Forward projection took ", time.time() - tic, " seconds,"
        else:
            # Ax
            b_hat = A.matvec(x) / b_weights 

        # DEBUGGING
        # if i == 1:
        #     b_debug = np.reshape(b_hat, (db.nt*db.nv, db.ns*db.nu))
        #     save_image("lf_" + str(i) + ".tif", b_debug);

        # STEP 2: Compute error between computed and observed light field images
        last_error = error
        if stabilization == "anscombe":
#            error = -b / anscombe( b_hat ) + 2
            error = b - b_hat
        else:
            error = b - b_hat

        # Get average of previous estimates
        if multiscale_smoothing:
            last_estimate_avg = multiscale_coefficient_mean(x, transform_type=transform_type)
        else:
            last_estimate_avg = np.mean(x)

        # Collect the unweighted error in light field space.
        if save_errors:
            iteration_error.append( error )

        # Add AMP adjustment.
        reweighted_error = error + (1.0/delta)*last_estimate_avg*last_error
        print "\t--> Last estimage average", last_estimate_avg
    
        # STEP 3: Back-project error onto the volume (or wavelet coefficients)
        tic = time.time()
        error_backprojected = A.rmatvec(reweighted_error) 
        error_backprojected /= x_weights
        print "\t--> Backward projection took ", time.time() - tic, " seconds,"

        if multiscale_smoothing:
            tic = time.time()

            # These are the multiscale coefficient errors. There is some reshaping that 
            # must go on to pass the volumes between lflib (column-major (y,x,z)) and 
            # R (row-major (x,y,z) if an R method is used. 
            if transform_type == "pyramidal_median":
                error_backprojected, vol_inds = multiscale_transform_3D(error_backprojected, 
                                                                        vol_shape=vol_shape,                                                                         
                                                                  wavelet_type=wavelet_type, 
                                                              transform_type=transform_type) 
            else:
                error_backprojected = error_backprojected.reshape(vol_shape).transpose((1,0,2))
                error_backprojected = error_backprojected.flatten(order='f')
                error_backprojected, vol_inds = multiscale_transform_3D(error_backprojected, 
                                                                        vol_shape=vol_shape,                                                                         
                                                                  wavelet_type=wavelet_type, 
                                                              transform_type=transform_type) 

            print "\t--> Multiscale transform took ", time.time() - tic, " seconds,"

        # Apply the update. 
        if multiscale_smoothing:
            # Update multiscale coefficients.
            update_norm = 0
            x, update_norm  = multiscale_coefficient_update(x, error_backprojected, alpha,
                                                            transform_type=transform_type)

            # Threshold updated coefs to impose nonnegativity and sparsity.
            tic = time.time()
            print "\t--> Thresholding multiscale coefficients..."
            threshold = [0.0, 0.0, 0.0, 0.0, 1e10, 1e10] 
#            threshold = sparsity
            x = multiscale_threshold(x, threshold=threshold, transform_type=transform_type)
            print "\t--> Multiscale thresholding took ", time.time() - tic, " seconds."
        else:
            # standard AMP update
            x_update = error_backprojected
            x += alpha * x_update

            # threshold element-wise to impose nonnegative and L1
            print "\t--> Sparsity:", sparsity
            x[x<sparsity]=0.0

        # CHECK FOR CONVERGENCE
        #
        # normalize MSE using input LF
        nrays = db.ns*db.nt*db.nu*db.nv
        residual_norm = np.linalg.norm(error) / nrays

        # normalize MSE using input LF
        nvoxels = db.nx*db.ny*db.nz
        if multiscale_smoothing:
            update_norm /= nvoxels
        else:
            update_norm = np.linalg.norm(alpha * x_update) / nvoxels

        toc = time.time()
        print ''
        print '\t[ AMP Iteration %d   (%0.2f seconds) ] ' % (i, toc-iteration_tic)
        print '\t      Residual Norm: %0.4g' % (residual_norm)
        print '\t        Update Norm: %0.4g               (tol = %0.2e)  ' % (update_norm, convergence_threshold)

        # check if convergence criteria met
        if i > 0 and residual_norm < convergence_threshold:
            break 

    if multiscale_smoothing:
        # Final projection back to volume space
        print "\t--> Final thresholding..."
        threshold = [0.0, 0.0, 0.0, 1e10, 1e10, 1e10] 
        x = multiscale_threshold(x, threshold=threshold, transform_type=transform_type, 
                                suppress_scales=[0])

        print "\t--> Preparing multiscale coefs for visualization..."
        multiscale_coefs_R = output_multiscale_coefs(x,J=J)
        multiscale_coefs = []
        for i in xrange(J):
            if standardize:
                scale = lf_max*np.asarray(multiscale_coefs_R[i])
            else:
                scale = np.asarray(multiscale_coefs_R[i])
            scale = scale.astype(np.float32)
            multiscale_coefs.append(scale)

        print "\t--> Final back projection to volume..."
        if transform_type == "pyramidal_median":
            x = np.asarray(inverse_multiscale_transform_3D(x,vol_inds,
                                                   transform_type=transform_type)).flatten()
        else:
            x = np.asarray(inverse_multiscale_transform_3D(x,vol_inds,
                                transform_type=transform_type)).transpose((1,0,2)).flatten()
#            x[x<0.0]=0.0

    x = inverse_anscombe(x)
    if standardize:
        vol = np.reshape(x*lf_max, (db.ny, db.nx, db.nz))
    else:
        vol = np.reshape(x, (db.ny, db.nx, db.nz))

    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    vol[0:db.supersample_factor, :, :] = 0.0
    vol[-db.supersample_factor:, :, :] = 0.0
    vol[:, 0:db.supersample_factor, :] = 0.0
    vol[:, -db.supersample_factor:, :] = 0.0    

    if multiscale_smoothing:
        print "Returning volume and multiscale coefficents..."
        return vol.astype(np.float32), multiscale_coefs
    else:
        print "Returning volume..."
        return vol.astype(np.float32), None


if __name__ == "__main__":
    pass

