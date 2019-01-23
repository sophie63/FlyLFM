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
from lflib.linear_operators import LightFieldOperator, RegularizedNormalEquationLightFieldOperator

#------------------------------------------------------------------------------------

class CgIterationFunctor:
    def __init__(self, algorithm_name, linear_operator, b, nvoxels, nrays, db):
        self.algorithm_name = algorithm_name
        self.iterations = 0
        self.prev_x = None
        self.last_time = 0
        self.linear_operator = linear_operator
        self.b = b
        self.nvoxels = nvoxels
        self.nrays = nrays
        self.db = db

    def iter_callback(self, x):
        toc = time.time()

        if self.prev_x != None:
            print '\t--> [ CG Iteration %d   (%0.2f seconds) ] ' % (self.iterations,
                                                                   toc - self.last_time)
            if self.linear_operator != None:
                b_hat = self.linear_operator.matvec(x)
                residual = self.b - b_hat
                residual_norm = np.linalg.norm(residual) / self.nrays
                print '\t      Residual Norm: %0.4g' % (residual_norm)
                update_norm = np.linalg.norm(self.prev_x-x) / self.nvoxels
                print '\t        Update Norm: %0.4g' % (update_norm)
                    
        
        self.last_time = toc
        self.prev_x = np.copy(x)
        self.iterations += 1

# ----------------------------------------------------------------------------------------
#                            CONJUGATE GRADIENT SOLVER
# ----------------------------------------------------------------------------------------

def conjugate_gradient_reconstruction(lfcal, lightfield, 
                                      convergence_threshold, max_iterations, regularization_lambda, 
                                      disable_gpu = False, gpu_id = 0):

    if lfcal.psf_db != None:
        db = lfcal.psf_db
    else:
        db = lfcal.rayspread_db

    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db, disable_gpu = disable_gpu, gpu_id = gpu_id)
    lfproj.set_premultiplier(lfcal.radiometric_correction)

    # Trim out entries in the light field that we are ignoring because
    # they are too close to the edge of the NA of the lenslet.  This
    # make our residuals more accurate below.
    lightfield = lfcal.mask_trimmed_angles(lightfield)

    # Form the b matrix 
    im_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)
    b = np.reshape(im_subaperture, (im_subaperture.shape[0]*im_subaperture.shape[1]))

    # Uncomment to swap between zero initial volume and focal stack
    # volume.  (Zero volume start seems to converge faster, though).
    vol_vec = np.zeros((db.ny*db.nx*db.nz), dtype=np.float32)
    #vol_vec = np.reshape(lfproj.backproject(lightfield), (db.nx * db.ny * db.nz))

    # Conjugate gradient requires a square A matrix, so we solve the
    # normal equations below, rather than the original problem Ax = b.
    from lflib.linear_operators import LightFieldOperator
    from lflib.linear_operators import NormalEquationLightFieldOperator
    from lflib.linear_operators import RegularizedNormalEquationLightFieldOperator

    # This reweighting factor should account for the non-stationary
    # variance of the shot-noise dominated data.
    vol_ones = db.ones_volume()
    lf_ones = lfproj.project(vol_ones)
    lf_ones_vec = np.reshape(lf_ones.asimage(LightField.TILED_SUBAPERTURE), (db.nrays))
    weight_idxs = np.nonzero(b != 0)
    weighting = np.zeros_like(b)
    weighting[weight_idxs] = 1.0 / (b[weight_idxs])
    weighting[np.nonzero(lf_ones_vec == 0.0)] = 0;
    
    A_operator = LightFieldOperator(lfproj, db)
    AtA_operator = NormalEquationLightFieldOperator(lfproj, db, weighting)
    A = LinearOperator((db.nt*db.nv*db.ns*db.nu, db.nx*db.ny*db.nz),
                       matvec=A_operator.matvec, rmatvec=A_operator, dtype='float')
    AtA = LinearOperator((db.nx*db.ny*db.nz, db.nx*db.ny*db.nz),
                         matvec=AtA_operator.matvec, rmatvec=AtA_operator.matvec, dtype='float')
    At_b = A_operator.rmatvec(weighting * b)  # b / b = 1

    print 'Calling poisson-weighted Conjugate Gradient solver...'
    from scipy.sparse.linalg import cg

    nrays = db.ns*db.nt*db.nu*db.nv
    nvoxels = db.nx*db.ny*db.nz
    iter_func = CgIterationFunctor("CG", A, b, nvoxels, nrays, db)
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

#---------------------------------------------------------------------------------------

if __name__ == "__main__":
    pass
