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

// Code for numerically integrating over the unit disk
#include "gauss_2d_sphere.h"

// The GNU Scientific Library (GSL) is used for sine integral function gsl_sf_Si()
#include <gsl/gsl_sf_expint.h>

// -----------------------        UTILITY FUNCTIONS       ----------------------------

inline double sinc(double x) { return sin(x)/x; }
inline double  unit_box(double x) { return (fabs(x) < 0.5) ? 1.0 : 0.0; }


//----------------------------------------------------------------------------------
//                     PROJECTIVE GEOMETRY HELPER ROUTINES
//----------------------------------------------------------------------------------

struct ProjectivePlane3D {
  float _coeffs[4];
  ProjectivePlane3D(float a, float b, float c, float x0, float y0, float z0) {
    _coeffs[0] = a; _coeffs[1] = b; _coeffs[2] = c; _coeffs[3] = -a*x0-b*y0-c*z0;
  }

  float operator[](int idx) const { return _coeffs[idx]; }
};

struct ProjectivePoint3D {
  float _coeffs[4];
  ProjectivePoint3D(float a, float b, float c, float d) {
    _coeffs[0] = a; _coeffs[1] = b; _coeffs[2] = c; _coeffs[3] = d;
  }

  float operator[](int idx) const { return _coeffs[idx]; }
};


// Compute the intersection of three projective planes in 3-space.
// 
// This closed-form solution was derived using Mathematica.
//
ProjectivePoint3D plane_intersect(ProjectivePlane3D const& p1, 
                                ProjectivePlane3D const& p2, 
                                ProjectivePlane3D const& p3) {

  float a1 = p1[0]; float b1 = p1[1]; float c1 = p1[2]; float d1 = p1[3];
  float a2 = p2[0]; float b2 = p2[1]; float c2 = p2[2]; float d2 = p2[3];
  float a3 = p3[0]; float b3 = p3[1]; float c3 = p3[2]; float d3 = p3[3];

  return ProjectivePoint3D ( (b3*c2*d1 - b2*c3*d1 - b3*c1*d2 + b1*c3*d2 + b2*c1*d3 - b1*c2*d3) /
                             (-a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3),
                             
                             (a3*c2*d1 - a2*c3*d1 - a3*c1*d2 + a1*c3*d2 + a2*c1*d3 - a1*c2*d3) /
                             (a3*b2*c1 - a2*b3*c1 - a3*b1*c2 + a1*b3*c2 + a2*b1*c3 - a1*b2*c3),
                             
                             (a3*b2*d1 - a2*b3*d1 - a3*b1*d2 + a1*b3*d2 + a2*b1*d3 - a1*b2*d3) /
                             (-a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3),
                      
                             1.0 );
}


//----------------------------------------------------------------------------------
//                            PIXEL TO ANGLE CONVERSION
//----------------------------------------------------------------------------------

//  Convert from discrete ray angle direction to an angle.  The
//  lenslet class implements this use the Abbe sine correction.
//  Returns a tuple (theta, phi) of angles relative to the optical
//  axis of the lenslet.
bool uv_to_angle(float &theta, float& phi, float u, float v, int nu, int nv, float na, 
                 float medium_index,  bool sine_correction = true) {

  // Still todo:
  // - Check (nu-1) for r_na below

  float r_na = (nu-1)/2.0;
  float r = sqrt( powf( (float(u)-r_na), 2.0) + powf( (float(v)-r_na), 2.0) );

  // If the ray is outside of the numerical aperture, we clamp it to
  // the numerical aperture and return false to indicate that we did
  // the clamping (just in case we want to throw away this ray, rather than clamping).
  bool returnval;
  if (r > r_na) {

    theta = -asin(na / medium_index);
    returnval = false;

  } else {

    double mag = 20;
    double pixel_size = 4.55e-6;
    double ulens_focal_length = 2433e-6;
    double medium_index = 1.33;
    double alpha_i = atan2(r * pixel_size, ulens_focal_length );

    if (sine_correction) {
      // theta = -asin(r / r_na * na / medium_index);  // Old Levoy / Oldenbourg style abbe sine correction
      theta = -asin(mag / medium_index * tan(alpha_i));
    } else {
      // theta = -r / r_na * asin(na / medium_index); // Old uncorrected
      theta = -atan(mag / medium_index * tan(alpha_i));
    }
    returnval = true;
  }

  // Compute phi, and handle the case where r = 0 (and phi therefore
  // does not matter, and tend to be computed as NaN...)
  phi = atan2( (float(v)-r_na), (float(u)-r_na) ) + M_PI;
  if (r == 0) 
    phi = 0;

  return returnval;
}

//----------------------------------------------------------------------------------
//                            ANTIALIASING KERNELS 
//----------------------------------------------------------------------------------


// -------------------------------------------
// Anti-aliasing for CIRCULAR lenslet profiles
// -------------------------------------------

struct CostFunctionMetadata {
  float m, n;  // Current offset of the filtering kernel
  int supersample_factor;  
};

const int NUM_LANCZOS_TAPS = 1;
static inline int cost_function_buffer() { return NUM_LANCZOS_TAPS; }

double lanczos_cost_fn(double x, double y, void* data) { 

  // Parse out metadata
  CostFunctionMetadata* meta = static_cast<CostFunctionMetadata*>(data);
  int supersample_factor = meta->supersample_factor;
  float m = meta->m;  float n = meta->n;

  double rho = sqrt( (x-m)*(x-m)+(y-n)*(y-n) );
  return 1.0/(M_PI*0.25) *
    sinc(supersample_factor*M_PI*rho) *
    sinc(supersample_factor*M_PI*rho/float(NUM_LANCZOS_TAPS)) * 
    unit_box(supersample_factor*rho/(2.0*NUM_LANCZOS_TAPS));
}



// -----------------------------------------------------------
// Anti-aliasing for SQUARE, TRUNCATED SPHERE lenslet profiles
// -----------------------------------------------------------


// This class is used to pass a 2D rayspread projection (on a specific
// z-plane) into these anti-aliasing filters.
struct RayspreadRectangle {
  double ll_x, ll_y;
  double ur_x, ur_y;

  RayspreadRectangle( double lower_left_x,  double lower_left_y,
                      double upper_right_x, double upper_right_y):
    ll_x(lower_left_x), ll_y(lower_left_y), 
    ur_x(upper_right_x), ur_y(upper_right_y) {}
};

// Lanczos Windowed Sinc filter (alpha = 1)
struct SquareLensletLanczos1ReconstructionFilter {
  static inline double integrate(RayspreadRectangle const& rect) {
    
    float a = std::max(rect.ll_x, -1.0);
    float b = std::min(rect.ur_x,  1.0);
    float c = std::max(rect.ll_y, -1.0);
    float d = std::min(rect.ur_y,  1.0);

    if (a == 0)
      a += 1e-2;
    if (b == 0)
      b += 1e-2;
    if (c == 0)
      c += 1e-2;
    if (d == 0)
      d += 1e-2;

    // If either of these conditions are violated, then the rayspread
    // does not overlap with this filter.
    if (a >= b || c >= d)
      return 0.0;

    float h_integral = -a+b-b*cos(2*a*M_PI) + a*cos(2*b*M_PI) + 2*a*b*M_PI*(gsl_sf_Si(2 * M_PI * b) - 
                                                                            gsl_sf_Si(2 * M_PI * a));
    float v_integral = -c+d-d*cos(2*c*M_PI) + c*cos(2*d*M_PI) + 2*c*d*M_PI*(gsl_sf_Si(2 * M_PI * d) - 
                                                                            gsl_sf_Si(2 * M_PI * c));
    return 1.0 / (4*a*b*c*d*powf(M_PI, 4)) * h_integral * v_integral;
  }  

  // This static method helps us know how many samples to buffer our
  // rayspread kernels by for this AA filter.
  static inline int buffer() { return 1; }
};


// ----------------------- DEBUGGING CODE ----------------------------

// void test_lanczos(int supersample_factor, int _m, int _n) {
//   CostFunctionMetadata m;
//   m.supersample_factor = supersample_factor;
//   m.n = _n;
//   m.m = _m;

//   const int width = 1024;
//   const int height = 1024;
  
//   cv::Mat test_image(width, height, CV_32FC1);
  
//   for (int i = 0; i < width; ++i) {
//     float x = float(i)/(width-1) * 2.0 - 1.0;
//     for (int j = 0; j < height; ++j) {
//       float y = float(j)/(height-1) * 2.0 - 1.0;
//       test_image.at<float>(i,j) = lanczos_cost_fn(x,y,&m);
//       if (test_image.at<float>(i,j) < 0) {
//         std::cout << x << " " << y << " " << test_image.at<float>(i,j) << "\n";
//       }
//     }
//   }
//   double minval, maxval;
//   cv::minMaxLoc(test_image, &minval, &maxval);
//   std::cout << "MINMAX: " << minval << " " << maxval << "\n";
//   test_image = test_image.mul(1.0/maxval * 255);
//   cv::minMaxLoc(test_image, &minval, &maxval);
//   std::cout << "MINMAX: " << minval << " " << maxval << "\n";
//   std::ostringstream ostr;
//   ostr << "test-" << supersample_factor << "-" << _m << "-" << _n << ".png";
//   std::cout << "Saving " << ostr.str() << "\n";
//   cv::imwrite(ostr.str(), test_image);
// };


//----------------------------------------------------------------------------------
//                       RAYSPREAD MAIN DRIVER ROUTINE
//----------------------------------------------------------------------------------

// rayspread()
//
// A few notes on coordinates systems:
// 
// Ray spreads are computed in a normalized coordinate system where
// each lenslet has a diameter of 1.0.  The center of the lenslet is
// (0,0).
// 
// 
float* rayspread(float z_depth, int u, int v, float nu, float nv, float z_spacing,
                 float pitch, float magnification, float na, 
                 float medium_index, bool sine_correction, 
                 int supersample_factor, float* radiometry,
                 int& kernel_h_offset, int& kernel_v_offset,
                 int& kernel_width, int& kernel_height,
		 int num_raysamples, int num_zsamples,
		 float pixel_fill_factor, float lenslet_fill_factor, 
                 int circular_ulens_profile) {

  // DEBUGGING
  // std::cout << "Testing lanczos3...\n";
  // for (int i = 0; i < supersample_factor; ++i) {
  //   for (int j = 0; j < supersample_factor; ++j) {
  //     test_lanczos(supersample_factor, i, j);
  //   }
  // }
  // exit(0);
  // /DEBUGGING

  const float PIXEL_FILL_FACTOR = pixel_fill_factor;
  const float LENSLET_FILL_FACTOR = lenslet_fill_factor;
  const int NUM_ZSAMPLES = num_zsamples;
  const int NUM_RAYSAMPLES = num_raysamples;
  
  // This is the size of the microlens in the sample plane
  float pitch_sampleplane = pitch / magnification;

  // Pixels behind the lenslet have finite size, so we must
  // average over "subpixel" u,v positions behind the lenslet to
  // get the full ray spread attributable to a single pixel.
  // This approximates each pixel by computing 10x10 subpixel
  // rayspreads between the range [ {u,v} - 0.45, {u,v} + 0.45 ]

  // ----------------------------------------------------------
  // FIRST PASS: 
  //
  // Determine the offsets and integer-aligned bounding box of the
  // kernel.  Here we only seek to the four "corners" of the pixel.
  // In the subsequent pass we'll sample the sub-rays of the pixel
  // with finer granularity.
  //
  float min_u = 1e99;  float min_v = 1e99;
  float max_u = -1e99; float max_v = -1e99;

  float zz = z_depth - z_spacing/2.0;
  for (float nzz = 0; nzz < NUM_ZSAMPLES; ++nzz) {
    float z_normalized = zz / pitch_sampleplane;

    float vv = v-PIXEL_FILL_FACTOR/2;
    for (int nvv = 0; nvv < 2; ++nvv) {    
      float uu = u-PIXEL_FILL_FACTOR/2;
      for (int nuu = 0; nuu < 2; ++nuu) {

        float theta, phi;
        bool within_na = uv_to_angle(theta, phi, uu, vv, nu, nv, na, medium_index, sine_correction);

        // Check if the result is within the NA of the objective.  If
        // not, we throw out this ray spread.
        if (!within_na) {
          kernel_width = 0;
          kernel_height = 0;
          kernel_h_offset = 0;
          kernel_v_offset = 0;
          return NULL;
        }
        
        ProjectivePoint3D p_center(float(z_normalized * tan(theta) * cos(phi)), 
                                   float(z_normalized * tan(theta) * sin(phi)),
                                   z_normalized, 1.0);
        ProjectivePoint3D p1(p_center[0]-LENSLET_FILL_FACTOR/2.0, 
                             p_center[1]-LENSLET_FILL_FACTOR/2.0, p_center[2], 1.0);
        ProjectivePoint3D p3(p_center[0]+LENSLET_FILL_FACTOR/2.0, 
                             p_center[1]+LENSLET_FILL_FACTOR/2.0, p_center[2], 1.0);
        //        std::cout << "a: " << p_center[0] << " " << p_center[1] << " " << p_center[2] << "\n";
            
        min_u = std::min(p1[0], min_u); 
        max_u = std::max(p3[0], max_u);
        min_v = std::min(p1[1], min_v);
        max_v = std::max(p3[1], max_v);
        
        uu += PIXEL_FILL_FACTOR;
      }
      vv += PIXEL_FILL_FACTOR;
    }
    zz += z_spacing / (NUM_ZSAMPLES-1);
  }

  // VERY IMPORTANT NOTE: It's important that these kernels are
  // aligned to the subsample_factor=1 grid.  The cl_backproject()
  // OpenCL kernel in volume.py expects this.
  //
  // TODO: Lanczos and other filters with a buffer() == 1 will violate
  // this!
  kernel_h_offset = floor(min_u)*supersample_factor - cost_function_buffer();
  kernel_v_offset = floor(min_v)*supersample_factor - cost_function_buffer();
  kernel_width    = ceil(max_u)*supersample_factor - floor(min_u)*supersample_factor + 
    2*cost_function_buffer() + 1;
  kernel_height   = ceil(max_v)*supersample_factor - floor(min_v)*supersample_factor + 
    2*cost_function_buffer() + 1;

  // Allocate memory for the kernel, and zero out that memory.
  float* kernel = new float[kernel_width * kernel_height];
  memset(kernel, 0.0, kernel_width*kernel_height*sizeof(float));

  // ----------------------------------------------------------
  // SECOND PASS: 
  //
  // Compute the kernel analytically, and then sample it (using the
  // appropirate raised cosine sampling filter to avoid
  // anti-aliasing).
  
  // Attenuation accounts for sampling of ray spread across 10
  // "sub-z-planes" and across NUM_RAYSAMPLES x NUM_RAYSAMPLES
  // "sub-angle" rays.
  float attenuation = NUM_RAYSAMPLES * NUM_RAYSAMPLES * NUM_ZSAMPLES; 

  // TODO: Add depth attenuation here (needs to be measured!)
  zz = z_depth - z_spacing/2.0;
  for (float nzz = 0; nzz < NUM_ZSAMPLES; ++nzz) {

    // Adjust the z-depth into "normalized" length units
    // where the lenslet is 1 unit x 1 unit in size.
    float z_normalized = zz / pitch_sampleplane;

    float vv = v-PIXEL_FILL_FACTOR/2;
    for (int nvv = 0; nvv < NUM_RAYSAMPLES; ++nvv) {
      float uu = u-PIXEL_FILL_FACTOR/2;
      for (int nuu = 0; nuu < NUM_RAYSAMPLES; ++nuu) {

        // RADIOMETRY
        //
        // Use the supplied nu x nv pixel radiometry (with linear
        // interpolation) to compute the radiometric correction for
        // this ray spread.
        //
        int u_0, u_1, v_0, v_1;
        if (vv - v <= 0) { v_0 = v  ; v_1 = v-1; } else { v_0 = v  ; v_1 = v+1; }
        if (uu - u <= 0) { u_0 = u  ; u_1 = u-1; } else { u_0 = u  ; u_1 = u+1; }

        // Clamp to make sure we don't index off the edge of the radiometry image.
        v_0 = std::min(std::max(0,v_0),int(nv)-1);  v_1 = std::min(std::max(0,v_1), int(nv)-1);
        u_0 = std::min(std::max(0,u_0),int(nu)-1);  u_1 = std::min(std::max(0,u_1), int(nu)-1);

        float r_00 = radiometry[ v_0  * int(nu) + u_0  ];
        float r_01 = radiometry[ v_0  * int(nu) + u_1 ];
        float r_10 = radiometry[ v_1 * int(nu) + u_0  ];
        float r_11 = radiometry[ v_1 * int(nu) + u_1 ];
        
        float dv = fabs(vv - v);  float du = fabs(uu - u);
        float radiometry_scale_factor = r_00*(1-dv)*(1-du) + r_10*dv*(1-du) + r_01*(1-dv)*du + r_11*du*dv;

        // Apply attenuation due to the fact that we are approximating
        // the overall rayspread by the sum of many discrete ray
        // spread parts.
        radiometry_scale_factor /= attenuation;

        // For dubugging: comment out above and uncomment this line to
        // turn off linear interp for radiometry.
        //  float radiometry_scale_factor = radiometry[ v * int(nu) +  u ];


        // PROJECTION & FILTERING
        //
        // Compute the center of the projected lenslet image, and then
        // sample it using the reconstruction filter.
        //
        float theta, phi;
        uv_to_angle(theta, phi, uu, vv, nu, nv, na, medium_index, sine_correction);
        ProjectivePoint3D p_center(float(z_normalized * tan(theta) * cos(phi)), 
                                   float(z_normalized * tan(theta) * sin(phi)),
                                   z_normalized, 1.0);

        for (int n = 0; n < kernel_height; ++n) {
          for (int m = 0; m < kernel_width; ++m) {

            // Multiply the rayspread corners by the
            // supersample_factor, and shift the rayspread so that it
            // is centered on the current voxel sample.
            float h_shift = float(kernel_h_offset + m + cost_function_buffer() ) / supersample_factor;
            float v_shift = float(kernel_v_offset + n + cost_function_buffer() ) / supersample_factor;
            
            // Compute the numerical integral of the reconstruction filter on
            // the unit disk using gaussian quadrature code by Pavel
            // Holoborodko (more info here: http://www.holoborodko.com/pavel/?page_id=1879)
            CostFunctionMetadata meta;
            meta.supersample_factor = supersample_factor;
            meta.m = h_shift;
            meta.n = v_shift; 

            if (circular_ulens_profile) {
              kernel[n * kernel_width + m] += radiometry_scale_factor * 
                gauss_product_2D_sphere(2 * supersample_factor,      // Quadrature order = 2x supersample factor
                                        lanczos_cost_fn, &meta,      // 
                                        LENSLET_FILL_FACTOR/2,       // Rayspread radius
                                        -p_center[0], -p_center[1]); // Rayspread center
            } else {
              h_shift = float(kernel_h_offset+m);
              v_shift = float(kernel_v_offset+n);
              ProjectivePoint3D p1(p_center[0]-LENSLET_FILL_FACTOR/2.0, 
                                   p_center[1]-LENSLET_FILL_FACTOR/2.0, p_center[2], 1.0);
              ProjectivePoint3D p3(p_center[0]+LENSLET_FILL_FACTOR/2.0, 
                                   p_center[1]+LENSLET_FILL_FACTOR/2.0, p_center[2], 1.0);
              RayspreadRectangle ray_rect(p1[0]*supersample_factor - h_shift, 
                                          p1[1]*supersample_factor - v_shift, 
                                          p3[0]*supersample_factor - h_shift, 
                                          p3[1]*supersample_factor - v_shift);
              kernel[n * kernel_width + m] += radiometry_scale_factor * 
                                        SquareLensletLanczos1ReconstructionFilter::integrate(ray_rect);
            }
          }
        }

        uu += PIXEL_FILL_FACTOR / (NUM_RAYSAMPLES-1);
      }
      vv += PIXEL_FILL_FACTOR / (NUM_RAYSAMPLES-1);
    }
    zz += z_spacing / (NUM_ZSAMPLES-1);
  }  

  // Adjust the kernel offsets to account for the buffer.  The
  // project() and backproject() methods will expect these offsets to
  // point to the unbuffered pixel location.
  kernel_h_offset += cost_function_buffer();
  kernel_v_offset += cost_function_buffer();
  return kernel;
}
