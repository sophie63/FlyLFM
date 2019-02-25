// __BEGIN_LICENSE__
//
// Copyright (C) 2010-2012 Stanford University.
// All rights reserved.
//
// __END_LICENSE__

#include <iostream>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <list>

#include <arroyo.h>

// -------------------------------------------------------------------------------------------

// light_field_psf()
//
// Uses the Arroyo library to produce a "splat" for a point at
// coordinates x,y,z in the volume.  The resulting intensity image
// takes into account diffraction from the objective and lenslet
// array.
// 
float* light_field_psf(float x, float y, float z, float um_per_slice,
                       int ns, int nt, int nu, int nv,
                       float ulens_focal_length, float ulens_pitch,
                       float objective_magnification, float objective_na, 
                       float medium_index, float wavelength, int supersample_factor, 
                       float pixel_fill_factor, float lenslet_fill_factor,
                       int& psf_width, int& psf_height) {

  // Half angles (computed from NA of objective and ulens array)
  double ulens_na = objective_magnification * atan2(ulens_pitch,(2*ulens_focal_length));
  double ulens_theta = asin(ulens_na / medium_index);
  double objective_theta = asin(objective_na / medium_index);
  
  // Effective optical parameters for the lenslet array (when
  // conjugated to the object side of the objective).
  double effective_pitch = ulens_pitch / objective_magnification * 1e-6;
  double effective_focal_length = (effective_pitch/2.0) / tan(ulens_theta);

  long final_wf_pix_per_lenslet = nu;
  long final_wf_pix_per_xform = nu;
  long num_lenslets = ns;

  const float pixel_scale = effective_pitch / final_wf_pix_per_lenslet;
  const float sim_size = pixel_scale * final_wf_pix_per_lenslet * num_lenslets;

  // If z = 0, much of the code below will fail.  We can "patch" this
  // up by checking if z is zero, and setting it ot a slightly
  // non-zero value if it is.
  if (z == 0) 
    z = 0.01e-6;

  // Create a plane wave
  Arroyo::three_frame wavefront_frame;
  Arroyo::diffractive_wavefront_header<double> dwfh(std::vector<long>(2,int(sim_size / 
                                                                            pixel_scale)),   // dimensions
                                                    wavefront_frame,                         // direction
                                                    wavelength, pixel_scale);

  // Create a square lenslet array
  std::vector<long> lenslet_axes(2,num_lenslets);
  Arroyo::square_lenslet_array sq_lnslt_array(lenslet_axes,
                                              effective_focal_length,
                                              effective_pitch,
                                              final_wf_pix_per_lenslet,
                                              final_wf_pix_per_xform);

  Arroyo::diffractive_wavefront<double> dwf(dwfh);

  // We spread an equal amount of light intensity over the area of the
  // aperture.
  double aperture_radius = fabs( z*tan(objective_theta) );
  // Baseline intensity for aperture radius equal to 0.5 * lenslet_pitch
  double I1 = 1.0/(M_PI*(effective_pitch/2.0)*(effective_pitch/2.0));  
  double k = 2*aperture_radius / effective_pitch;
  // Adjusted intensity for aperture with radius k * baseline_radius
  double I2 = I1 / k;
  
  //  dwf += std::complex<double>(1.0/(M_PI*aperture_radius*aperture_radius), 0.0);
  dwf += std::complex<double>(I2, 0.0);

  // Set the curvature of the wavefront to the reciprocal of the distance to the point source.
  dwf.set_wavefront_curvature(1.0/z);

  // Create an initial aperture of the appropriate size to trim out
  // light outside the NA of the objective.  
  Arroyo::circular_aperture circ_ap( 2*aperture_radius );
  circ_ap.transform(dwf);

  // TODO: Add radiometry here.

  // Propagate the wavefront through the microlens array.  By
  // default, arroyo will propagate the light to the lenslet array's
  // focal plane.
  Arroyo::three_translation trans(Arroyo::three_vector(x, y, 0.0, wavefront_frame));
  trans.transform(sq_lnslt_array);
  sq_lnslt_array.transform(dwf);

  // Extract the data
  psf_width = final_wf_pix_per_lenslet*ns;
  psf_height = final_wf_pix_per_lenslet*nt;
  float* psf = new float[psf_width * psf_height];
  memset(psf, 0.0, psf_width*psf_height*sizeof(float));

  std::vector<long> axes = dwf.get_axes();
  for (int j = 0; j < axes[1]; ++j) {
    for (int i = 0; i < axes[0]; ++i) {
      int n = i*axes[0]+j;
      std::complex<double> d = dwf.data(n);
      psf[j * psf_width + i] = real(d * conj(d));
    }
  }
  return psf;
}

