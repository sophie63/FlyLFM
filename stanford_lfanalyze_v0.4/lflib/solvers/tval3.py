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

#-------------------------------------------------------------------------------
#                               TVAL3 SOLVER
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------

if __name__ == "__main__":
    pass
