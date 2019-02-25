# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
cimport numpy as np
np.import_array()

from libc.stdint cimport *
from libcpp cimport bool
from cython.operator cimport dereference as deref

# This code is needed to transparently convert from python strings to
# stl strings.
cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()

cdef extern from "stdlib.h":
    void free(void* ptr)

#---------------------------------------------------------------------------------------
#                                  C++ WRAPPERS
#---------------------------------------------------------------------------------------

cdef extern from "rayspread.h":
    float* rayspread(float z, int u, int v, float nu, float nv,
                     float z_spacing, 
                     float pitch, float magnification, float na,
                     float medium_index, bool sine_correction,
                     int supersampling, float* radiometry,
                     int& kernel_h_offset, int& kernel_v_offset,
                     int& kernel_width, int& kernel_height,
 		     int num_raysamples, int num_zsamples,
		     float pixel_fill_factor, float lenslet_fill_factor,
                     bool circular_ulens_profile)


cdef extern from "longobject.h":
    float* long_object_compensation( int nu, int nv, int ns, int nt, 
                                 int num_slices, float um_per_slice,
                                 float z_center,
                                 float magnification,
                                 float pitch,
                                 float medium_index,
                                 float na )

# cdef extern from "physical_model.h":
#     float* light_field_psf(float x, float y, float z, float um_per_slice,
#                            int ns, int nt, int nu, int nv,
#                            float ulens_focal_length, float ulens_pitch,
#                            float objective_magnification, float objective_na, 
#                            float refractive_index, float wavelength, int supersample_factor, 
#                            float pixel_fill_factor, float lenslet_fill_factor,
#                            int& psf_width, int& psf_height)


#---------------------------------------------------------------------------------------
#                                PYTHON CLASSES
#---------------------------------------------------------------------------------------

def compute_rayspread(lenslet_array, u, v, z, z_spacing,
                      np.ndarray[np.float32_t,ndim=2] radiometry,
                      sine_correction = True, supersample_factor = 1,
		      num_raysamples = 10, num_zsamples = 5, pixel_fill_factor = 0.9,
                      ulens_fill_factor = 0.7, ulens_profile = 'rect'):
    cdef int kernel_width, kernel_height, kernel_h_offset, kernel_v_offset
    cdef float* kernel
    if ulens_profile == 'rect':
        circular_ulens_profile = False
    elif ulens_profile == 'circ':
        circular_ulens_profile = True
    else:
        print 'ERROR: \'',ulens_profile,'\' is an unrecognized microlens profile.  Choose \'rect\' or \'circ\''

    kernel = rayspread(z, u, v, lenslet_array.nu, lenslet_array.nv, z_spacing,
                       lenslet_array.ulens_pitch, lenslet_array.objective_magnification,
                       lenslet_array.objective_na,
                       lenslet_array.medium_index, sine_correction, supersample_factor,
                       <float*> radiometry.data,
                       kernel_h_offset, kernel_v_offset, kernel_width, kernel_height,
		       num_raysamples, num_zsamples, pixel_fill_factor, ulens_fill_factor,
                       int(circular_ulens_profile))

    # Handle the case where we are returned an empty (i.e. nonexistent) rayspread.
    if kernel_width == 0 or kernel_height == 0:
        return (np.zeros((0,0)), kernel_h_offset, kernel_v_offset)

    cdef np.npy_intp intp_size = kernel_width * kernel_height
    cdef np.ndarray newarr = np.PyArray_SimpleNewFromData(1, &intp_size, np.NPY_FLOAT, <void *>kernel)
    result = newarr.copy()
    free(kernel)
    return (result.reshape((kernel_height, kernel_width)), kernel_h_offset, kernel_v_offset)

# def compute_light_field_psf(x, y, z, z_spacing, supersample_factor, num_lenslets_in_psf,
#                             lenslet_array, np.ndarray[np.float32_t,ndim=2] radiometry):
#     cdef int psf_width, psf_height
#     cdef float* psf
#     psf = light_field_psf(x, y, z, z_spacing,
#                           num_lenslets_in_psf, num_lenslets_in_psf, lenslet_array.nu, lenslet_array.nv,
#                           lenslet_array.ulens_focal_length, lenslet_array.ulens_pitch,
#                           lenslet_array.objective_magnification, lenslet_array.objective_na,
#                           lenslet_array.medium_index, 525e-9, # Wavelength set to green for now...
#                           supersample_factor, lenslet_array.ulens_fill_factor, lenslet_array.pixel_fill_factor,
#                           psf_width, psf_height)

#     # Handle the case where we are returned an empty (i.e. nonexistent) rayspread.
#     if psf_width == 0 or psf_height == 0:
#         return np.zeros((0,0))

#     cdef np.npy_intp intp_size = psf_width * psf_height
#     cdef np.ndarray newarr = np.PyArray_SimpleNewFromData(1, &intp_size, np.NPY_FLOAT, <void *>psf)
#     result = newarr.copy()
#     free(psf)
#     return result.reshape(psf_height, psf_width)

def compute_long_object_compensation( lfcal ):
    raydb = lfcal.rayspread_db
        
    nu = raydb.nu
    nv = raydb.nv
    ns = raydb.ns
    nt = raydb.nt
        
    num_slices    = lfcal.num_z_slices
    um_per_slice  = lfcal.um_per_z_slice
    z_center      = lfcal.z_center
    magnification = lfcal.magnification
    pitch         = lfcal.array_pitch
    medium_index  = lfcal.medium_index
    na            = lfcal.na
		
    cdef float* data
    data = long_object_compensation( nu, nv, ns, nt, num_slices, um_per_slice,
		                             z_center, magnification, pitch, medium_index, na )
    
    cdef np.npy_intp intp_size = nv * nu * ns * nt
    cdef np.ndarray newarr = np.PyArray_SimpleNewFromData(1, &intp_size, np.NPY_FLOAT, <void *>data)
    result = newarr.copy()
    free( data )
    return result.reshape(nv*nt, nu*ns)
