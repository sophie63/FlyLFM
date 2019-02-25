// __BEGIN_LICENSE__
// Zahid Hossain
// Copyright (C) 2010-2012 Stanford University.
// All rights reserved.
//
// __END_LICENSE__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#include "rayspread.h"
#include "longobject.h"

/***
    Zahid Hossain
    Assumption: This function assumes that the angles are not close to 90 degrees
    which should be the case most of the time.
*/
float long_object_ratio_one_side( float tan_angle, float min_h, float max_h, float min_z, float max_z)
{   
        
    float h_max = max_z * tan_angle; 
    float h_min = min_z * tan_angle;

	if( ((h_max < min_h) && (h_min < min_h )) ||
        ((h_max > max_h) && (h_min > max_h )) )
        return 0;
    
    if( (h_max >= min_h) && (h_max <= max_h ) &&
        (h_min >= min_h) && (h_min <= max_h ))
        return 1;
    
    
    float height = max_z - min_z;
    float hh,zz;
    float cot_angle = 1.0 / tan_angle;
    
    if( (h_max >= min_h) && (h_max <= max_h ) ){
        hh = std::min(std::max( h_min, min_h ),max_h);        
        zz = hh * cot_angle;
        return (max_z - zz) / height;
    }
    else{
        hh = std::min(std::max( h_max, min_h ),max_h);        
        zz = hh * cot_angle;
        return (zz - min_z) / height;
    }
}
        

float long_object_ratio( float tan_theta, float tan_phi, int s, int t, int ns, int nt, 
                         float x_range[], float y_range[], float z_range[], float pitch )
{       
    float x = (s - (ns-1)*0.5) * pitch;
    float y = (t - (nt-1)*0.5) * pitch;
    
    float xx_range[2] = {x_range[0] - x, x_range[1] - x};
    float yy_range[2] = {y_range[0] - y, y_range[1] - y};
            
    float phi_ratio   = long_object_ratio_one_side( tan_phi, yy_range[0], yy_range[1], z_range[0], z_range[1] );
    float theta_ratio = long_object_ratio_one_side( tan_theta, xx_range[0], xx_range[1], z_range[0], z_range[1] );    
    
        //Maybe we need to compute this master_ratio later
    float master_ratio = 1.0;
    float ratio = std::min(theta_ratio,phi_ratio) * master_ratio;   
    return ratio*ratio*ratio;
}

/*****
   Zahid Hossain
   This computes a mask, having the same dimension as the light-field image
   organized as TILED_SUBAPERTURE type, compensate for the long-object problem. 
   
   We are using a very simple approximation based on the following paper:
   
   Wei Xua, Fang Xua, Mel Jones, Bettina Keszthelyi, John Sedat, David Agard, Klaus Mueller, 
   "High-performance iterative electron tomography reconstruction with long-object
   compensation using graphics processing units (GPUs)," Journal of structural biology, 2010
   
   Except we made a very minor modification of raising the compensation ratio to the power of 3
   because our ray-spreads are shaped like pyramids, so has a growing volume component in it.
      
*/
float* long_object_compensation( int nu, int nv, int ns, int nt, 
                                 int num_slices, float um_per_slice,
                                 float z_center,
                                 float magnification,
                                 float pitch,
                                 float medium_index,
                                 float na )
{
        
#define LF_PIXEL(u,v,s,t) ( (v)*nu*ns*nt + (u)*ns + (t) * ns * nu + (s) )
#define LENSLET_PIXEL(u,v) ( (v)*nu + u )
    
    float *data = (float*) malloc ( nu * nv * ns * nt * sizeof(float) );
    float actual_pitch = pitch / magnification;
    
    float x_range[2] = { - (ns-1) * 0.5 * actual_pitch , (ns-1) * 0.5 * actual_pitch };
    float y_range[2] = { - (nt-1) * 0.5 * actual_pitch , (nt-1) * 0.5 * actual_pitch };
    float z_range[2] = { z_center - (num_slices - 1) * um_per_slice * 0.5, z_center + (num_slices - 1) * um_per_slice * 0.5 };
    
    float* theta_map = (float*) malloc( nu * nv * sizeof(float) );
    float* phi_map   = (float*) malloc( nu * nv * sizeof(float) );
    
    // Precomputing the tangent of angles ( For the angles it's always using Abbe Sine correction )
    for( int v = 0; v < nv; ++v ){
        for( int u = 0; u < nu; ++u ){
          float theta, phi;
          uv_to_angle( theta, phi,u,v,nu,nv,na,medium_index,true );
          theta_map[ LENSLET_PIXEL(u,v) ] = tan(theta);
          phi_map  [ LENSLET_PIXEL(u,v) ] = tan(phi);
        }
    }
        
        // Loop is ordered such that pixels are updated one scanline at a time
        // for cache efficiency     
    for( int v = 0; v < nv; ++v ){         
        for( int t = 0; t < nt; ++t ){
            for( int u = 0; u < nu; ++u ){
                for( int s = 0; s < ns; ++s ){
                    data [ LF_PIXEL(u,v,s,t) ] = long_object_ratio ( theta_map[ LENSLET_PIXEL(u,v) ], 
                                                                             phi_map[ LENSLET_PIXEL(u,v) ], 
                                                                             s, t, ns, nt, 
                                                                             x_range, y_range, z_range,
                                                                             actual_pitch );
                }
            }
        }
    } 
    
    free( theta_map );
    free( phi_map );
    
    return data;
}


