// __BEGIN_LICENSE__
//
// Copyright (C) 2010-2012 Stanford University.
// All rights reserved.
//
// __END_LICENSE__

float* light_field_psf(float x, float y, float z, float um_per_slice,
                       int ns, int nt, int nu, int nv,
                       float ulens_focal_length, float ulens_pitch,
                       float objective_magnification, float objective_na, 
                       float refractive_index, float wavelength, int supersample_factor, 
                       float pixel_fill_factor, float lenslet_fill_factor,
                       int& psf_width, int& psf_height);

