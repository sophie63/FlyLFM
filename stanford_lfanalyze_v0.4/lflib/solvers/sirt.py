# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2013 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
import time

from lflib.imageio import save_image

# ----------------------------------------------------------------------------------------
#                             SIRT ITERATIVE SOLVER
# ----------------------------------------------------------------------------------------

def sirt_reconstruction(A, b, alpha, tol, max_iter):
    """
    Deconvolve the volume using the Simultaneous Iterative
    Reconstruction Techinque (SIRT), adapted for light field
    microscopy.For background reading on SIRT and related algorithms,
    we recommend the following references:
    
    A. C. Kak and Malcolm Slaney, "Principles of Computerized Tomographic Imaging,"
    Society of Industrial and Applied Mathematics, 2001.
    
    Jiang, Ming, and Ge Wang. "Development of iterative algorithms for image reconstruction."
    Journal of XRay Science and Technology 10.2 (2002) : 77-86.
    """

    # The A operator represents a large, sparse matrix that has dimensions [ nrays x nvoxels ]
    nrays = A.shape[0]
    nvoxels = A.shape[1]

    # The Neo has 2e- RMS read noise, but I measured it to be 2.81
    # DN so we use that here.
    # readnoise_variance = 7.84;  
    # test = np.maximum(b - 600, 1.0)
    preconditioner = 1.0 /np.sqrt(b + 1.0) 
    #    A_operator.diagonalizer = preconditioner
    #    A = A_operator.as_linear_operator(nrays, nvoxels)
    b_orig = b.copy()
    b *= preconditioner

    # Create the SIRT weight volume and light field.  These are
    # created by projecting a volume containing all ones (to create
    # the weight lightfield), and then back-projecting a light field
    # with all ones (to create the weight volume).
    #
    # x_weights = A.T * b_ones;    b_weights = A * x_ones
    import time

    tic = time.time()
    b_weights = preconditioner * A.matvec(np.ones((nvoxels), dtype=np.float32))
    print '\t--> Time for one forward projection: %0.2f seconds.' % (time.time()-tic)

    tic = time.time()
    x_weights = A.rmatvec(preconditioner * np.ones((nrays), dtype=np.float32))
    print '\t    Time for one back projection: %0.2f seconds.' % (time.time()-tic)

    # Make sure that any zero valued weights are set to nonzero so
    # that they don't lead to division by zero below. We then
    # normalize the starting volume using the volume weights.
    min_bweight = b_weights[np.nonzero(b_weights != 0.0)].min()
    min_xweight = x_weights[np.nonzero(x_weights != 0.0)].min()
    b_weights[np.nonzero(b_weights < min_bweight)] = min_bweight;
    x_weights[np.nonzero(x_weights < min_xweight)] = min_xweight;

    # Start with a zero'd out initial volume.
    x = np.zeros(nvoxels, dtype=np.float32)

    # --------------------------------------------------------------------
    
    iteration_residuals = []
    for i in range(max_iter):
        tic = time.time()

        # In each iteration, forward and backproject residual from all views at once, then update the volume
        #
        # STEP 1: forward projection of volume to create sub-aperture images.
        Ax = A.matvec(x)
        b_hat = preconditioner * Ax

        # STEP 2: Compute residual between computed and observed sub-aperture images
        residual = b - b_hat

        # if (i == 9):
        #     ns = 93; nt = 77; nu = 28; nv = 28;
        #     save_image('debug_residual_weighted.tif', np.reshape(residual, (nt*nv, ns*nu)), dtype=np.float32)
        #     save_image('debug_residual_raw.tif', np.reshape(b_orig-Ax, (nt*nv, ns*nu)), dtype=np.float32)
        #     save_image('debug_bhat.tif', np.reshape(Ax, (nt*nv, ns*nu)), dtype=np.float32)
        #     save_image('debug_b.tif', np.reshape(b_orig, (nt*nv, ns*nu)), dtype=np.float32)
            
        iteration_residuals.append( residual )

        # the reweighted residual
        reweighted_residual = residual / b_weights

        # STEP 3: back-project residual onto the volume and update
        residual_backprojected = A.rmatvec(preconditioner * reweighted_residual)
        x_update = (residual_backprojected / x_weights)

        # Apply the update
        x += alpha * x_update

        # Enforce non-negativity constraint
        x[x<0]=0

        # CHECK FOR CONVERGENCE
        #
        # normalize MSE using input LF
        residual_norm = np.linalg.norm(residual) / nrays

        # normalize MSE using input LF
        update_norm = np.linalg.norm(alpha * x_update) / nvoxels

        toc = time.time()
        print '\t--> [ SIRT Iteration %d   (%0.2f seconds) ] ' % (i, toc-tic)
        print '\t      Residual Norm: %0.4g' % (residual_norm)
        print '\t        Update Norm: %0.4g               (tol = %0.2e)  ' % (update_norm, tol)

        # check if convergence criteria met
        if i > 0 and update_norm < tol:
            break 

    return (x.astype(np.float32), iteration_residuals)

#-----------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#EOF
