# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np

#---------------------------------------------------------------------------------------
#                                  UTILITY METHODS
#---------------------------------------------------------------------------------------

def subaperture2lenslet(lf_subaperture, nu, nv, ns, nt):
    """
    Turn the light field representation inside-out so that each
    lenslet image is stored in adjacent memory.  This function returns
    an nu*ns x nv*nt image where ns x nt lenslet images are adjacent
    to each other in the image.
    """
    lf_lenslet = np.zeros(lf_subaperture.shape, dtype = lf_subaperture.dtype)
    nCols = lf_lenslet.shape[1] / nu
    nRows = lf_lenslet.shape[0] / nv
    for v in range(nv):
        for u in range(nu):
            lf_lenslet[v:lf_lenslet.shape[0]:nv,
                       u:lf_lenslet.shape[1]:nu] = lf_subaperture[v*nRows:v*nRows+nRows,
                                                                  u*nCols:u*nCols+nCols]
    return lf_lenslet

def lenslet2subaperture(lf_lenslet, nu, nv, ns, nt):
    """
    Turn the light field representation inside-out so that each
    sub-aperture image is stored in adjacent memory.  This function
    returns an nu*ns x nv*nt image where nu x nv aperture images are
    adjacent to each other in the image.
    """
    lf_subaperture = np.zeros(lf_lenslet.shape, dtype = lf_lenslet.dtype)
    nCols = lf_lenslet.shape[1] / nu
    nRows = lf_lenslet.shape[0] / nv
    for v in range(nv):
        for u in range(nu):
            lf_subaperture[v*nRows:v*nRows+nRows,
                           u*nCols:u*nCols+nCols] = lf_lenslet[v:lf_lenslet.shape[0]:nv,
                                                               u:lf_lenslet.shape[1]:nu]
    return lf_subaperture


def lenslet_image_coords_to_uvst(vec, nu, nv):
    # Converts a series of image coordinates in a TILED_LENSLET light
    # field image stored in the rows of a vector to light field
    # coordinates.  IMPORTANT NOTE: [m, n] are expected to be in numpy
    # order (i.e. [row, col])
    m = vec[:,0]
    n = vec[:,1]

    u = np.mod(n, nu);
    v = np.mod(m, nv);
    s = np.floor(np.divide(n,nu))
    t = np.floor(np.divide(m,nv))

    return np.vstack((u, v, s, t)).T.astype(np.int32)

def uvst_to_lenslet_image_coords(vec, nu, nv):
    # Converts a series of light field coords stored in the rows of a
    # vector to image pixel coordinates in a TILED_LENSLET light field
    # image.  IMPORTANT NOTE: [m, n] are returned in numpy order
    # (i.e. [row, col])
    u = vec[:,0]
    v = vec[:,1]
    s = vec[:,2]
    t = vec[:,3]

    return np.vstack((t*nv + v, s*nu + u)).T.astype(np.int32)
    

#---------------------------------------------------------------------------------------
#                                  LIGHTFIELD CLASS
#---------------------------------------------------------------------------------------

class LightField(object):
    """
    The lightfield class stores a light field with u x v ray angles
    and s x t pixel positions.
    """

    # Enumeration: represents how raw lightfield data is stored in
    # memory.
    TILED_SUBAPERTURE = 0
    TILED_LENSLET     = 1

    # Constructor.
    def __init__(self, data, nu, nv, ns, nt, representation = TILED_LENSLET):
        self._data = data
        self._representation = representation

        self.nu = nu; self.nv = nv
        self.ns = ns; self.nt = nt


    # Return raw lightfield data as an image.
    def asimage(self, representation = TILED_LENSLET):
        """Returns the raw data for the light field in the requested representation."""
        if representation == self._representation:
            return self._data
        elif self._representation == self.TILED_LENSLET:
            return lenslet2subaperture(self._data, self.nu, self.nv, self.ns, self.nt)
        elif self._representation == self.TILED_SUBAPERTURE:
            return subaperture2lenslet(self._data, self.nu, self.nv, self.ns, self.nt)

    def subaperture(self, u, v):
        if self._representation == self.TILED_SUBAPERTURE:
            return self._data[self.nt*v:self.nt*v+self.nt, self.ns*u:self.ns*u+self.ns]
        elif self._representation == self.TILED_LENSLET:
            return self._data[v::self.nv, u::self.nu]
            
