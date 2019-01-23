// __BEGIN_LICENSE__
//
// Copyright (C) 2010-2012 Stanford University.
// All rights reserved.
//
// __END_LICENSE__

#include <utility>

bool uv_to_angle(float &theta, float &phi, float u, float v, 
                 int nu, int nv, float na, 
                 float medium_index, bool sine_correction = true);
    
float* rayspread(float z, int u, int v, float nu, float nv, float z_spacing,
                 float pitch, float magnification, float na, 
                 float medium_index, bool sine_correction,
                 int supersampling, float* radiometry,
                 int& kernel_h_offset, int& kernel_v_offset,
                 int& kernel_width, int& kernel_height,
		 int num_raysamples, int num_zsamples,
		 float pixel_fill_factor, float lenslet_fill_factor, 
                 int circular_ulens_profile);

