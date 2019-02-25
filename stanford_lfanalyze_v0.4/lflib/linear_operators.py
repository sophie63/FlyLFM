# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2013 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
import time

from lflib.lightfield import LightField
from scipy.sparse.linalg.interface import LinearOperator

#------------------------------------------------------------------------------------
#                 LINEAR OPERATORS FOR LIGHT FIELD RECONSTRUCTION
#------------------------------------------------------------------------------------

class ScaledLinearOperator(object):
    def __init__(self, A, scale_factor):
        self.A = A
        self.scale_factor = scale_factor

    def matvec(self, x ):
        return self.A.matvec(x) * self.scale_factor

    def rmatvec(self, x ):
        return self.A.rmatvec(x) * self.scale_factor

#------------------------------------------------------------------------------------

class LightFieldOperator(object):
    def __init__(self, sirt, db):
        self.sirt = sirt
        self.db = db
        self.diagonalizer = None

    def matvec(self, vol_vec ):
        vol = np.reshape(vol_vec, (self.db.ny, self.db.nx, self.db.nz))
        im = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        b = np.reshape(im, (im.shape[0]*im.shape[1]))
        if self.diagonalizer == None:
            return b
        else:
            return self.diagonalizer * b

    def rmatvec(self, lf_vec ):
        if self.diagonalizer == None:
            lf = np.reshape(lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        else:
            lf = np.reshape(self.diagonalizer*lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                          representation = LightField.TILED_SUBAPERTURE))
        return np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))

    def as_linear_operator(self, nrays, nvoxels):
        return LinearOperator((nrays, nvoxels),
                              matvec=self.matvec,
                              rmatvec=self.rmatvec,
                              dtype='float')

#------------------------------------------------------------------------------------
class LightFieldOperatorWithIntercept(object):
    def __init__(self, sirt, db, intercept):
        self.sirt = sirt
        self.db = db
        self.intercept = intercept
        self.diagonalizer = None

    def matvec(self, vol_vec ):
        volume_term = vol_vec[0:-1]
        backrgound_term = vol_vec[-1]
        print 'BACKGROUND TERM: ', backrgound_term

        vol = np.reshape(volume_term, (self.db.ny, self.db.nx, self.db.nz))
        im = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        b = np.reshape(im, (im.shape[0]*im.shape[1]))
        b += (backrgound_term * self.intercept)
        if self.diagonalizer == None:
            return b
        else:
            return self.diagonalizer * b

    def rmatvec(self, lf_vec ):
        if self.diagonalizer == None:
            lf = np.reshape(lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        else:
            lf = np.reshape(self.diagonalizer*lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                          representation = LightField.TILED_SUBAPERTURE))
        vol_vec =  np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))
        return np.append(vol_vec, np.dot(self.intercept.T, lf_vec))

    def as_linear_operator(self, nrays, nvoxels):
        return LinearOperator((nrays, nvoxels+1),
                              matvec=self.matvec,
                              rmatvec=self.rmatvec,
                              dtype='float')

class LightFieldOperatorWithFullIntercept(object):
    def __init__(self, sirt, db):
        self.sirt = sirt
        self.db = db
        self.diagonalizer = None

    def matvec(self, vol_vec ):
        volume_term = vol_vec[0:self.db.nvoxels]
        background_term = vol_vec[self.db.nvoxels:]
        # print 'BACKGROUND TERM: ', background_term.mean(), np.median(background_term)

        vol = np.reshape(volume_term, (self.db.ny, self.db.nx, self.db.nz))
        im = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        b = np.reshape(im, (im.shape[0]*im.shape[1]))
        b += background_term
        if self.diagonalizer == None:
            return b
        else:
            return self.diagonalizer * b

    def rmatvec(self, lf_vec ):
        if self.diagonalizer == None:
            lf = np.reshape(lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        else:
            lf = np.reshape(self.diagonalizer*lf_vec, (self.db.nt*self.db.nv, self.db.ns*self.db.nu))
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                          representation = LightField.TILED_SUBAPERTURE))
        vol_vec =  np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))
        return np.append(vol_vec, lf_vec)

    def as_linear_operator(self):
        return LinearOperator((self.db.nrays, self.db.nvoxels+self.db.nrays),
                              matvec=self.matvec,
                              rmatvec=self.rmatvec,
                              dtype='float')

#------------------------------------------------------------------------------------

class NormalEquationLightFieldOperator(object):
    def __init__(self, sirt, db, weighting_matrix = None):
        self.sirt = sirt
        self.db = db
        if weighting_matrix != None:
            self.weighting_matrix = np.reshape(weighting_matrix, (db.nv*db.nt, db.nu*db.ns))
        else:
            self.weighting_matrix = None

    def matvec(self, vol_vec ):
        vol = np.reshape(vol_vec.astype(np.float32), (self.db.ny, self.db.nx, self.db.nz))
        lf = self.sirt.project(vol).asimage(representation = LightField.TILED_SUBAPERTURE)
        if self.weighting_matrix != None:
            lf *= self.weighting_matrix
        vol = self.sirt.backproject(LightField(lf, self.db.nu, self.db.nv, self.db.ns, self.db.nt,
                                          representation = LightField.TILED_SUBAPERTURE))
        return np.reshape(vol, (self.db.nx * self.db.ny * self.db.nz))

#------------------------------------------------------------------------------------

class RegularizedNormalEquationLightFieldOperator(object):
    def __init__(self, sirt, db, regularization_lambda):
        self.sirt = sirt
        self.db = db
        self.regularization_lambda = regularization_lambda

    def matvec(self, vol_vec ):
        input_vol = np.reshape(vol_vec.astype(np.float32), (self.db.ny, self.db.nx, self.db.nz))
        lf = self.sirt.project(input_vol).asimage(representation = LightField.TILED_SUBAPERTURE)
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

#------------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#EOF
