# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import sys
import numpy as np
from lflib.lightfield import LightField
from lflib.imageio import save_image
from scipy.ndimage import filters

try:
    import pyopencl as cl
    import pyopencl.array as cl_array
    import pyopencl.characterize as cl_characterize
    LFLIB_HAVE_OPENCL = True
except ImportError:
    LFLIB_HAVE_OPENCL = False

# Utility function
extract = lambda x, y: dict(zip(x, map(y.get, x)))

# ------------------------------------------------------------------------------
#                               OPENCL KERNELS
# ------------------------------------------------------------------------------

src = '''

// Enable 64-bit double support
#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

__kernel void cl_backproject_rayspread(__read_only image2d_t lf,
                                 __global float* volume,
                                 int lf_rows, int lf_cols,
                                 __global float* psf,
                                 __constant int* filterRows,
                                 __constant int* filterCols,
                                 __constant int* filterRowOffset,
                                 __constant int* filterColOffset,
                                 int supersample_factor,
                                 sampler_t sampler) {
  
  // Store each work-item unique row and column
  int x = get_global_id(0);        int y = get_global_id(1); 
  int s = x / supersample_factor;  int t = y / supersample_factor;

  // Bail out if the kernel is out of bounds
  if (x >= K_NX || y >= K_NY)
    return;

  // This offset comes into play when supersample_factor > 1
  int x_bump = (x - s * supersample_factor);
  int y_bump = (y - t * supersample_factor);
  
  // Iterator for the filter
  int base_address = 0;
  
  for (int z_idx = 0; z_idx < K_NZ; ++z_idx) {
    float sum = 0.0f;
    int2 coords;

    // Iterate the filter rows
    for(int j = y_bump; j < filterRows[z_idx]; j+=supersample_factor) {
      coords.y = (y - j - filterRowOffset[z_idx])/supersample_factor;

      // Iterate over the filter columns
      int filterIdx = base_address + filterCols[z_idx]*j + x_bump;
      for(int i = x_bump; i < filterCols[z_idx]; i+=supersample_factor) {
        coords.x = (x - i - filterColOffset[z_idx])/supersample_factor;

        // Read a pixel from the image. A single channel image
        // stores the pixel in the 'x' coordinate of the returned vector.
        float4 pixel = read_imagef(lf, sampler, coords);
        sum += pixel.x * psf[filterIdx];
        filterIdx += supersample_factor;
      }
    }

    // Copy the data to the output image if the work-item is in bounds
    volume[y*K_NX*K_NZ + x*K_NZ + z_idx] += sum;
    base_address += (filterCols[z_idx] * filterRows[z_idx]);
  }
}

__kernel void cl_project_rayspread(__read_only image2d_t vol_slice,
                         __global float* subaperture,
                         __global float* psf,
                         __constant int* filterRows,
                         __constant int* filterCols,
                         __constant int* filterRowOffset,
                         __constant int* filterColOffset,
                         __constant int* u_coords,
                         __constant int* v_coords,
                         int num_rays,
                         int supersample_factor,
                         sampler_t sampler) {

  // Store each work-item unique row and column
  int s = get_global_id(0); int t = get_global_id(1);

  // Bail out if the kernel is out of bounds
  if (s >= K_NS || t >= K_NT)
    return;

  // Iterator for the filter
  int base_address = 0;

  for (int r = 0; r < num_rays; ++r) {
    float sum = 0.0f;
    int2 coords;

    // Iterate the filter rows
    for(int j = 0; j < filterRows[r]; ++j) {
      coords.y = t*supersample_factor + j + filterRowOffset[r];

      // Iterate over the filter columns
      int filterIdx = base_address + filterCols[r]*j;
      for(int i = 0; i < filterCols[r]; ++i) {
        coords.x = s*supersample_factor + i + filterColOffset[r];

        // Read a pixel from the image. A single channel image
        // stores the pixel in the 'x' coordinate of the returned vector.
        float4 pixel = read_imagef(vol_slice, sampler, coords);
        sum += pixel.x * psf[filterIdx++];
      }
    }

    // Copy the data to the output light field
    int u = u_coords[r];
    int v = v_coords[r];
    subaperture[v*K_NT*K_NU*K_NS + t*K_NU*K_NS + u*K_NS + s] += sum;
    base_address += (filterCols[r] * filterRows[r]);
  }
  
}


__kernel void cl_project_wavespread(__read_only image2d_t vol_slice,
                             __global double* subaperture_im,
                             __global int* u_coords,
                             __global int* v_coords,
                             __global int* s_coords,
                             __global int* t_coords,
                             __global float* coefficients,
                             int num_coefficients,
                             int supersample_factor,
                             int dx, int dy,
                             sampler_t sampler) {

  // Create local storage to speed up running summations.
  __private double sum_buf[K_NU * K_NV];
//  __private float sum_total[K_NU * K_NV];
//  __private float num_good[K_NU * K_NV];
//  __private float num_total[K_NU * K_NV];
  for (int v = 0; v < K_NV; ++v) {
    for (int u = 0; u < K_NU; ++u) {
      sum_buf[v*K_NU+u] = 0.0;
//      sum_total[v*K_NU+u] = 0.0;
//      num_good[v*K_NU+u] = 0.0;
//      num_total[v*K_NU+u] = 0.0;
    }
  }
  
  // Store each work-item unique row and column
  int s = get_global_id(0); int t = get_global_id(1);

  // Bail out if the kernel is out of bounds
  if (s >= K_NS || t >= K_NT)
    return;

  // Iterate over the psf coordinates and coefficients
  for (int i = 0; i < num_coefficients; ++i) {

    // Grab the appropriate pixel from the volume.  Here we assume
    // that a value of s_coords = 0 and t_coords = 0 refer to the
    // lenslet at the center of the volume.
    //
    int2 coords;
    coords.x = (s - s_coords[i]) * supersample_factor + dx;
    coords.y = (t - t_coords[i]) * supersample_factor + dy;

    float4 pixel = read_imagef(vol_slice, sampler, coords);

    // Copy the data to the output light field
    int u = u_coords[i];
    int v = v_coords[i];
    //num_total[v*K_NU+u] += 1;
    //if (coords.x >= 0 && coords.x < K_NS && coords.y >= 0 && coords.y < K_NT) {
    //num_good[v*K_NU+u] += 1;
    //}
    
    sum_buf[v*K_NU+u] += pixel.x * coefficients[i];
  }

  for (int v = 0; v < K_NV; ++v) {
    for (int u = 0; u < K_NU; ++u) {
      subaperture_im[v*K_NT*K_NU*K_NS + t*K_NU*K_NS + u*K_NS + s] += sum_buf[v*K_NU+u];
    }
  }
}

__kernel void cl_backproject_wavespread(__read_only image2d_t lf,
                                        __global double* volume,
                                        __global int* u_coords,
                                        __global int* v_coords,
                                        __global int* s_coords,
                                        __global int* t_coords,
                                        __global float* coefficients,
                                        int num_coefficients,
                                        int z, int dx, int dy,
                                        int supersample_factor,
                                        sampler_t sampler) {
  
  // Store each work-item unique row and column
  int s = get_global_id(0);
  int t = get_global_id(1);
  int x = s * supersample_factor + dx; 
  int y = t * supersample_factor + dy; 

  // Bail out if the kernel is out of bounds
  if (x >= K_NX || y >= K_NY || x < 0 || y < 0)
    return;

  double sum = 0.0;

  // Iterate over the psf coordinates and coefficients
  for (int i = 0; i < num_coefficients; ++i) {

    int2 coords;

    //float r = sqrt( (float) ((s - K_NS / 2) + (t - K_NT / 2)) );
    //float theta = atan2((float)(s - K_NS / 2), (float)(t - K_NT / 2));

    //float du = r * cos(theta);
    //float dv = r * sin(theta);

    // Check to make sure the s index is in bounds.
    int s_idx = s + s_coords[i];
    int t_idx = t + t_coords[i];
    if (s_idx < 0 || s_idx >= K_NS || t_idx < 0 || t_idx >= K_NT)
      continue;

    // Grab the appropriate pixel from the light field.  Here we assume
    // that a value of s_coords = 0 and t_coords = 0 refer to the
    // lenslet at the center of the volume.
    //
    coords.y = t_idx * K_NS + s_idx;
    coords.x = v_coords[i] * K_NU + u_coords[i];
    float4 pixel = read_imagef(lf, sampler, coords);

    // Copy the data to the output light field
    sum += pixel.x * coefficients[i];
  }
  volume[z*K_NX*K_NY + y*K_NX + x] += sum;
}

'''

# ------------------------------------------------------------------------------
#                        UTILITY CLASSES & FUNCTIONS
# ------------------------------------------------------------------------------

class OptimizedRayDatabase(object):
    '''
    A utility class that converts a ray database into a form where
    spreads for a given u,v are laid out adjacent in memory.  This
    makes them easy to ship off to the GPU for backprojection
    '''
    UV_ADJACENT = 1
    Z_ADJACENT = 2
    
    def __init__(self, rayspreads, nu, nv, nz, layout = UV_ADJACENT):
        self.layout = layout
        self.psfs = {}
        self.kernel_h_offsets = {}
        self.kernel_v_offsets = {}
        self.kernel_widths = {}
        self.kernel_heights = {}
        self.u_coords = {}
        self.v_coords = {}
        self.z_coords = {}

        if self.layout == self.UV_ADJACENT:
            for u in range(nu):
                for v in range(nv):
                    self.psfs[(u,v)] = []
                    self.kernel_h_offsets[(u,v)] = []
                    self.kernel_v_offsets[(u,v)] = []
                    self.kernel_widths[(u,v)] = []
                    self.kernel_heights[(u,v)] = []
                    self.u_coords[(u,v)] = []
                    self.v_coords[(u,v)] = []
                    self.z_coords[(u,v)] = []

                    # This performance shortcut avoids processing the rays
                    # that are outside of the NA of the objective.  These rays
                    # are vignetting and highly aberrated, and aren't so good
                    # to use in any case.
                    if (np.sqrt(np.power(u-(nu-1.0)/2.0,2) + np.power(v-(nv-1.0)/2.0,2)) >= nv/2.0):
                        continue

                    for z in range(nz):
                        try:
                            spread = rayspreads[(z,u,v)]
                            kernel_h_offset = spread[1]
                            kernel_v_offset = spread[2]
                            psf = spread[0]
                            kernel_width = psf.shape[1]
                            kernel_height = psf.shape[0]
                            self.add((u,v),u,v,z,psf,kernel_h_offset, kernel_v_offset, kernel_width, kernel_height)
                        except KeyError:
                            # Some keys are not in the ray database.
                            # If this is the case, we simply move onto
                            # the next key.
                            pass

                        
        elif self.layout == self.Z_ADJACENT:
            for z in range(nz):
                self.psfs[(z)] = []
                self.kernel_h_offsets[(z)] = []
                self.kernel_v_offsets[(z)] = []
                self.kernel_widths[(z)] = []
                self.kernel_heights[(z)] = []
                self.u_coords[(z)] = []
                self.v_coords[(z)] = []
                self.z_coords[(z)] = []

                for u in range(nu):
                    for v in range(nv):

                        # This performance shortcut avoids processing the rays
                        # that are outside of the NA of the objective.  These rays
                        # are vignetting and highly aberrated, and aren't so good
                        # to use in any case.
                        if (np.sqrt(np.power(u-(nu-1.0)/2.0,2) + np.power(v-(nv-1.0)/2.0,2)) >= nv/2.0):
                            continue

                        try:
                            spread = rayspreads[(z,u,v)]
                            kernel_h_offset = spread[1]
                            kernel_v_offset = spread[2]
                            psf = spread[0]
                            kernel_width = psf.shape[1]
                            kernel_height = psf.shape[0]
                            self.add((z),u,v,z,psf,kernel_h_offset, kernel_v_offset, kernel_width, kernel_height)
                        except KeyError:
                            # Some keys are not in the ray database.
                            # If this is the case, we simply move onto
                            # the next key.
                            pass

    def add(self, idx, u, v, z, psf, h_offset, v_offset, width, height):
        self.psfs[idx] += np.reshape(psf, (width*height)).tolist()
        self.kernel_h_offsets[idx].append(h_offset)
        self.kernel_v_offsets[idx].append(v_offset)
        self.kernel_widths[idx].append(width)
        self.kernel_heights[idx].append(height)
        self.u_coords[idx].append(u)
        self.v_coords[idx].append(v)
        self.z_coords[idx].append(z)

def roundUp(value, multiple):
    ''' Determine how far past the nearest multiple of the value.
    This is useful for padding image dimensions to nice powers of 16,
    which helps promote efficient memory access on the GPU.'''
    
    remainder = value % multiple;
    if remainder != 0:
        value += (multiple-remainder);
    return value;

def padVector(vec, multiple, val = 0.0):
    '''
    Pad a vector so that it has an integer multiple of 'multiple'
    entries.  This is useful when doing vector memory access on the
    GPU.
    '''
    newsize = roundUp(len(vec), multiple)
    padded_vec = np.zeros(newsize, dtype=vec.dtype)
    padded_vec[0:len(vec)] = vec
    return padded_vec

# ------------------------------------------------------------------------------
#                         LIGHT FIELD PROJECTION CLASS
# ------------------------------------------------------------------------------

class LightFieldProjection(object):

    def __init__(self, rayspread_db = None, psf_db = None, disable_gpu = False, gpu_id = None):

        self.rayspread_db = rayspread_db
        self.psf_db = psf_db
        if self.psf_db != None:
            db = self.psf_db
        else:
            db = self.rayspread_db

        # Premultiplier is optionaly applied after forward
        # projection, and before back projection.
        #
        # Postmultiplier is optionaly applied before forward
        # projection, and after back projection.
        #
        self.premultiplier = None
        self.postmultiplier = None
        
        if LFLIB_HAVE_OPENCL and not disable_gpu:
            # Set up OpenCL
            platform = cl.get_platforms()[0]
            from pyopencl.tools import get_gl_sharing_context_properties
            
            if sys.platform == "darwin":
                self.cl_ctx = cl.Context() 
            elif gpu_id is not None:
                self.cl_ctx = cl.Context(properties=[(cl.context_properties.PLATFORM, platform)],
                                         devices = [platform.get_devices()[gpu_id]])
            else:
                self.cl_ctx = cl.Context(properties=[(cl.context_properties.PLATFORM, platform)])

            # This poor-man's template substitution allows us to
            # hard-code some values in the kernel.
            preprocessed_src = src.replace('K_NU', str(db.nu))
            preprocessed_src = preprocessed_src.replace('K_NV', str(db.nv))
            preprocessed_src = preprocessed_src.replace('K_NS', str(db.ns))
            preprocessed_src = preprocessed_src.replace('K_NT', str(db.nt))
            preprocessed_src = preprocessed_src.replace('K_NX', str(db.nx))
            preprocessed_src = preprocessed_src.replace('K_NY', str(db.ny))
            preprocessed_src = preprocessed_src.replace('K_NZ', str(db.nz))

            self.cl_prog = cl.Program(self.cl_ctx, preprocessed_src).build()
            self.cl_queue = cl.CommandQueue(self.cl_ctx)

            # Set up OpenCL
            self.WORKGROUP_X = 16
            self.WORKGROUP_Y = 16
            self.WORKGROUP_Z = 1

            if self.psf_db != None:
                self.backproject = self.backproject_wavespread_gpu
                self.project = self.project_wavespread_gpu
            else:
                self.backproject = self.backproject_rayspread_gpu
                self.project = self.project_rayspread_gpu

        else:  # No OpenCL
            if self.psf_db != None:
                print "WARNING: your platform does not seem to support OpenCL.  There is no project/backproject implementations for wavespreads on the CPU.   Exiting."
                raise SystemExit
            else:
                print "WARNING: your platform does not seem to support OpenCL.  Using SIRT CPU implementation."
                self.backproject = self.backproject_cpu
                self.project = self.project_cpu

    def set_premultiplier(self, premultiplier):
        self.premultiplier = premultiplier

    def set_postmultiplier(self, postmultiplier):
        self.postmultiplier = postmultiplier

    # -----------------------------------------------------------------------------
    #                          PSF & SPLAT ROUTINES
    # -----------------------------------------------------------------------------

    def psf( self, uv = None ):
        """
        Use the rayspread database to build a stack illustrating the point
        spread function that is formed by combining all rays in the database.
        """

        rayspreads = self.rayspread_db.rayspreads
        z_coords = self.rayspread_db.z_coords

        nx = self.rayspread_db.nx
        ny = self.rayspread_db.ny
        nz = len(self.rayspread_db.z_coords)

        volume = np.zeros((ny, nx, nz), dtype='float32')
        cx = nx / 2
        cy = ny / 2

        for spread_key in rayspreads:
            spread = rayspreads[spread_key]
            ray_psf = spread[0]
            kernel_h_offset = spread[1]
            kernel_v_offset = spread[2]
            kernel_width = ray_psf.shape[1]
            kernel_height = ray_psf.shape[0]

            # Compute the z index.  
            z = spread_key[0]
            u = spread_key[1]
            v = spread_key[2]

            if (uv is None):
                volume[cy+kernel_v_offset:cy+kernel_v_offset+kernel_height,
                       cx+kernel_h_offset:cx+kernel_h_offset+kernel_width, z] += ray_psf
            elif (u == uv[0] and v == uv[1]):
                volume[cy+kernel_v_offset:cy+kernel_v_offset+kernel_height,
                       cx+kernel_h_offset:cx+kernel_h_offset+kernel_width, z] += ray_psf

        # Renormalize to the range [0.0, 1.0]
        volume /= volume.max()
        return volume

    def splat( self, delta_coords ):
        """
        Use the rayspread database to build a light field corresponding to a 
        delta function placed at the given coordinates by calling the 
        project function.

        Important note: this function expect you to provide a
        coordinate in the form:

        (row, col, z)  --  i.e. (y, x, z)
        """
        nx = self.rayspread_db.nx
        ny = self.rayspread_db.ny
        nz = len(self.rayspread_db.z_coords)

        volume = np.zeros((ny, nx, nz), dtype='float32')
        volume[delta_coords] = 1.0
 
        lf = self.project_cpu( volume, delta_coords[2] )
        return lf

    # ============================================================================================
    #                              RAYSPREAD PROJECT/BACKPROJECT
    # ============================================================================================

    def project_rayspread_gpu( self, volume, zslice = None ):
        """
        Using the given rayspread database, compute a light field by
        projecting the volume using rayspreads.
        
        Returns: A lightfield as a tiled set of sub-aperture
        (i.e. pinhole) images.
        """

        rayspreads = self.rayspread_db.rayspreads
        z_coords = self.rayspread_db.z_coords

        nu = self.rayspread_db.nu; nv = self.rayspread_db.nv;
        ns = self.rayspread_db.ns; nt = self.rayspread_db.nt;
        nx = self.rayspread_db.nx; ny = self.rayspread_db.ny; nz = self.rayspread_db.nz;

        # Create an empty lightfield image.
        lf = np.zeros((self.rayspread_db.nv * nt, self.rayspread_db.nu * ns), dtype = np.float32)

        # Upload volume slices to the GPU
        vol_slices = {}
        samp = cl.Sampler(self.cl_ctx, False, cl.addressing_mode.CLAMP, cl.filter_mode.NEAREST)
        for z in range(volume.shape[2]):
            vol_slice = volume[:,:,z].astype(np.float32)
            vol_slices[z] = cl.image_from_array(self.cl_ctx, vol_slice.copy(), 1, 'r')

        # Create a set of empty subaperture images to accumulate data into
        subaperture_buf = cl_array.zeros(self.cl_queue, (nt*nv, ns*nu), dtype=np.float32)

        # Set static kernel arguments
        kern = self.cl_prog.cl_project_rayspread
        kern.set_arg(1, subaperture_buf.data)
        kern.set_arg(10, np.uint32(self.rayspread_db.supersample_factor))
        kern.set_arg(11, samp)

        # Workgroup and global sizes
        localSize = (self.WORKGROUP_X, self.WORKGROUP_Y)
        globalSize = (roundUp(ns,self.WORKGROUP_X),
                      roundUp(nt,self.WORKGROUP_Y))

        # Concatenate all of the rayspreads into one long vector.
        optimized_raydb = OptimizedRayDatabase(rayspreads, self.rayspread_db.nu,
                                       self.rayspread_db.nv, self.rayspread_db.nz,
                                       layout = OptimizedRayDatabase.Z_ADJACENT)


        for z in range(nz):
            if len(optimized_raydb.kernel_widths[(z)]) == 0:
                continue

            psfs_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                 hostbuf=np.array(optimized_raydb.psfs[(z)], dtype=np.float32))

            kernel_h_offsets_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf=np.array(optimized_raydb.kernel_h_offsets[(z)], dtype=np.int32))

            kernel_v_offsets_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf=np.array(optimized_raydb.kernel_v_offsets[(z)], dtype=np.int32))

            kernel_widths_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                          hostbuf=np.array(optimized_raydb.kernel_widths[(z)], dtype=np.int32))

            kernel_heights_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                          hostbuf=np.array(optimized_raydb.kernel_heights[(z)], dtype=np.int32))

            u_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                     hostbuf=np.array(optimized_raydb.u_coords[(z)], dtype=np.int32))

            v_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                     hostbuf=np.array(optimized_raydb.v_coords[(z)], dtype=np.int32))

            kern.set_arg(0, vol_slices[z])
            kern.set_arg(2, psfs_buf)
            kern.set_arg(3, kernel_heights_buf) 
            kern.set_arg(4, kernel_widths_buf)
            kern.set_arg(5, kernel_v_offsets_buf)
            kern.set_arg(6, kernel_h_offsets_buf)
            kern.set_arg(7, u_coords_buf)
            kern.set_arg(8, v_coords_buf)
            kern.set_arg(9, np.uint32(len(optimized_raydb.kernel_widths[(z)]))) # num_rays

            # Execute the kernel
            cl.enqueue_nd_range_kernel(self.cl_queue, kern, globalSize, localSize)

        # Download the result, and place the data for this
        # subaperture into the lightfield.
        lf = subaperture_buf.get(self.cl_queue)

        self.cl_queue.finish()
        return LightField(lf, self.rayspread_db.nu, self.rayspread_db.nv,
                          self.rayspread_db.ns, self.rayspread_db.nt,
                          representation = LightField.TILED_SUBAPERTURE)

    def backproject_rayspread_gpu( self, light_field):
        """
        Using the given rayspread database, compute focused images at the
        z_depths included in the supplied rayspread_db.
        
        Returns: A volume as a 3D numpy matrix that is indexed in [y,x,z]
        order.
        """
        rayspreads = self.rayspread_db.rayspreads

        nu = self.rayspread_db.nu
        nv = self.rayspread_db.nv
        nx = self.rayspread_db.nx
        ny = self.rayspread_db.ny
        nz = self.rayspread_db.nz
        
        # Create an empty volume
        volume_buf = cl_array.zeros(self.cl_queue, (ny,nx,nz), dtype=np.float32)

        # Grab the light field image in sub-aperture mode
        lf_textures = {}
        samp = cl.Sampler(self.cl_ctx, False, cl.addressing_mode.CLAMP, cl.filter_mode.NEAREST)
        for u in range(self.rayspread_db.nu):
            for v in range(self.rayspread_db.nv):
                lf = light_field.subaperture(u,v).astype(np.float32)
                lf_textures[(u,v)] = cl.image_from_array(self.cl_ctx, lf, 1, 'r')

        # Set static kernel arguments
        kern = self.cl_prog.cl_backproject_rayspread
        kern.set_arg(1, volume_buf.data)
        kern.set_arg(2, np.uint32(lf.shape[0]))
        kern.set_arg(3, np.uint32(lf.shape[1]))
        kern.set_arg(9, np.uint32(self.rayspread_db.supersample_factor))
        kern.set_arg(10, samp)

        localSize = (self.WORKGROUP_X, self.WORKGROUP_Y)
        globalSize = (roundUp(nx,self.WORKGROUP_X),
                      roundUp(ny,self.WORKGROUP_Y))

        # Concatenate all of the rayspreads into one long vector.
        uv_raydb = OptimizedRayDatabase(rayspreads, self.rayspread_db.nu,
                                        self.rayspread_db.nv, self.rayspread_db.nz,
                                        layout = OptimizedRayDatabase.UV_ADJACENT)

        for u in range(nu):
            for v in range(nv):
                if len(uv_raydb.kernel_widths[(u,v)]) == 0:
                    continue

                psfs_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                     hostbuf=np.array(uv_raydb.psfs[(u,v)], dtype=np.float32))

                kernel_h_offsets_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                                 hostbuf=np.array(uv_raydb.kernel_h_offsets[(u,v)], dtype=np.int32))

                kernel_v_offsets_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                                 hostbuf=np.array(uv_raydb.kernel_v_offsets[(u,v)], dtype=np.int32))

                kernel_widths_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                              hostbuf=np.array(uv_raydb.kernel_widths[(u,v)], dtype=np.int32))

                kernel_heights_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                              hostbuf=np.array(uv_raydb.kernel_heights[(u,v)], dtype=np.int32))

                kern.set_arg(0, lf_textures[(u,v)])
                kern.set_arg(4, psfs_buf)
                kern.set_arg(5, kernel_heights_buf)
                kern.set_arg(6, kernel_widths_buf)
                kern.set_arg(7, kernel_v_offsets_buf)
                kern.set_arg(8, kernel_h_offsets_buf)

                # Execute the kernel
                cl.enqueue_nd_range_kernel(self.cl_queue, kern, globalSize, localSize)
                
        # Download the result, and place the data for this subaperture into the lightfield
        volume = volume_buf.get(self.cl_queue)
        self.cl_queue.finish()
        return np.reshape(volume, (ny, nx, nz))

    # ============================================================================================
    #                              WAVESPREAD PROJECT/BACKPROJECT
    # ============================================================================================

    def project_wavespread_gpu( self, volume, zslice = None ):
        """
        Using the given psf database, compute a light field by
        projecting the volume onto the light field sensor.
        
        Returns: A lightfield as a tiled set of sub-aperture
        (i.e. pinhole) images.
        """

        if self.postmultiplier != None:
            vol = volume * self.postmultiplier
        else:
            vol = volume

        psf_coordinates = self.psf_db.psf_coordinates
        psf_coefficients = self.psf_db.psf_coefficients
        z_coords = self.psf_db.z_coords

        nu = self.psf_db.nu; nv = self.psf_db.nv;
        ns = self.psf_db.ns; nt = self.psf_db.nt;
        nx = self.psf_db.nx; ny = self.psf_db.ny; nz = self.psf_db.nz;
        supersample_factor = self.psf_db.supersample_factor

        # Create an empty lightfield image.
        lf = np.zeros((nv * nt, nu * ns), dtype = np.float64)

        # Upload volume slices to the GPU
        vol_slices = {}
        samp = cl.Sampler(self.cl_ctx, False, cl.addressing_mode.CLAMP, cl.filter_mode.NEAREST)
        for z in range(vol.shape[2]):
            vol_slice = vol[:,:,z].astype(np.float32)
            vol_slices[z] = cl.image_from_array(self.cl_ctx, vol_slice.copy(), 1, 'r')
 
        # Create a set of empty subaperture images to accumulate data into
        subaperture_buf = cl_array.zeros(self.cl_queue, (nt*nv, ns*nu), dtype=np.float64)

        # Set static kernel arguments
        kern = self.cl_prog.cl_project_wavespread
        kern.set_arg(1, subaperture_buf.data)
        kern.set_arg(8, np.uint32(self.psf_db.supersample_factor))
        kern.set_arg(11, samp)

        # Workgroup and global sizes
        localSize = (self.WORKGROUP_X, self.WORKGROUP_Y)
        globalSize = (roundUp(ns, self.WORKGROUP_X),
                      roundUp(nt, self.WORKGROUP_Y))

        for z in range(nz):
            for x in range(supersample_factor):
                for y in range(supersample_factor):
                    coords = psf_coordinates[(x,y,z)]
                    coefficients = psf_coefficients[(x,y,z)]

                    # Extract the sorted coordinates. Copy here ensure
                    # "single segment" memory buffer, which seems to
                    # be required by cl.Buffer() below.
                    u_coords = coords[:,0].copy()
                    v_coords = coords[:,1].copy()
                    s_coords = coords[:,2].copy()
                    t_coords = coords[:,3].copy()

                    u_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = u_coords)
                    v_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = v_coords)
                    s_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = s_coords)
                    t_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = t_coords)
                    coefficients_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                                 hostbuf=coefficients)

                    kern.set_arg(0, vol_slices[z])
                    kern.set_arg(2, u_coords_buf)
                    kern.set_arg(3, v_coords_buf)
                    kern.set_arg(4, s_coords_buf)
                    kern.set_arg(5, t_coords_buf)
                    kern.set_arg(6, coefficients_buf)
                    kern.set_arg(7, np.uint32(coefficients.shape[0])) # num_coefficients
                    kern.set_arg(9, np.uint32(x))
                    kern.set_arg(10, np.uint32(y))

                    # Execute the kernel
                    cl.enqueue_nd_range_kernel(self.cl_queue, kern, globalSize, localSize)

        # Download the result, and place the data for this
        # subaperture into the lightfield.
        self.cl_queue.finish()
        lf_im = subaperture_buf.get(self.cl_queue)
        lf = LightField(lf_im, self.psf_db.nu, self.psf_db.nv,
                        self.psf_db.ns, self.psf_db.nt,
                        representation = LightField.TILED_SUBAPERTURE)

        if self.premultiplier != None:
            light_field_im = lf.asimage(LightField.TILED_LENSLET) * self.premultiplier
            lf = LightField(light_field_im, nu, nv, ns, nt,
                            representation = LightField.TILED_LENSLET)

        self.cl_queue.finish()
        return lf

    def backproject_wavespread_gpu( self, light_field):
        """
        Using the given psf database, compute focused images at the
        z_depths included in the supplied psf_db.
        
        Returns: A volume as a 3D numpy matrix that is indexed in [y,x,z]
        order.
        """
        psf_coordinates = self.psf_db.psf_coordinates
        psf_coefficients = self.psf_db.psf_coefficients
        
        z_coords = self.psf_db.z_coords

        nu = self.psf_db.nu
        nv = self.psf_db.nv
        ns = self.psf_db.ns
        nt = self.psf_db.nt
        nx = self.psf_db.nx
        ny = self.psf_db.ny
        nz = self.psf_db.nz
        supersample_factor = self.psf_db.supersample_factor

        # Optionally apply radiometry correction
        if self.premultiplier != None:
            light_field_im = light_field.asimage(LightField.TILED_LENSLET) * self.premultiplier
            light_field = LightField(light_field_im, nu, nv, ns, nt,
                                     representation = LightField.TILED_LENSLET)

        # Create an empty volume
        volume_buf = cl_array.zeros(self.cl_queue,
                                    (ny*nx*nz),
                                    dtype=np.float64)

        # Grab the light field image in sub-aperture mode
        samp = cl.Sampler(self.cl_ctx, False, cl.addressing_mode.CLAMP, cl.filter_mode.NEAREST)

        # Store the PSF in a CL 2d image texture, which will give us
        # some speedup in the algorithm since it caches texture accesses.
        lf_texture_volume = np.zeros((nt*ns, nu*nv), dtype=np.float32)
        for u in range(nu):
            for v in range(nv):
                lf_texture_volume[:,v*nu+u] = np.reshape(light_field.subaperture(u,v), (ns*nt), order='C')
        lf_texture_2d = cl.image_from_array(self.cl_ctx, lf_texture_volume.copy(), 1, 'r')

        # Set static kernel arguments
        kern = self.cl_prog.cl_backproject_wavespread
        kern.set_arg(0, lf_texture_2d)
        kern.set_arg(1, volume_buf.data)
        kern.set_arg(11, np.uint32(supersample_factor))
        kern.set_arg(12, samp)

        localSize = (self.WORKGROUP_X, self.WORKGROUP_Y)
        globalSize = (roundUp(nx/supersample_factor,self.WORKGROUP_X),
                      roundUp(ny/supersample_factor,self.WORKGROUP_Y))

        for z in range(nz):
            for x in range(supersample_factor):
                for y in range(supersample_factor):
                    coords = psf_coordinates[(x,y,z)]
                    coefficients = psf_coefficients[(x,y,z)]

                    # Extract the sorted coordinates. Copy here ensure
                    # "single segment" memory buffer, which seems to
                    # be required by cl.Buffer() below.
                    u_coords = coords[:,0].copy()
                    v_coords = coords[:,1].copy()
                    s_coords = coords[:,2].copy()
                    t_coords = coords[:,3].copy()
                    
                    u_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = u_coords)
                    v_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = v_coords)
                    s_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = s_coords)
                    t_coords_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                             hostbuf = t_coords)
                    coefficients_buf = cl.Buffer(self.cl_ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                                                 hostbuf=coefficients)

                    kern.set_arg(2, u_coords_buf)
                    kern.set_arg(3, v_coords_buf)
                    kern.set_arg(4, s_coords_buf)
                    kern.set_arg(5, t_coords_buf)
                    kern.set_arg(6, coefficients_buf)
                    kern.set_arg(7, np.uint32(coefficients.shape[0])) # num_coefficients
                    kern.set_arg(8, np.uint32(z))  
                    kern.set_arg(9, np.uint32(x))  # dx
                    kern.set_arg(10, np.uint32(y))  # dy

                    # Execute the kernel
                    cl.enqueue_nd_range_kernel(self.cl_queue, kern, globalSize, localSize)
                
        # Download the result, and place the data for this subaperture into the lightfield
        volume = volume_buf.get(self.cl_queue)
        self.cl_queue.finish()
        vol = np.transpose(np.reshape(volume, (nz, ny, nx)),  (1, 2, 0))
        
        if self.postmultiplier != None:
            return vol * self.postmultiplier
        else:
            return vol



    # ----------------------------------------------------------------------------
    #                          CPU IMPLEMENTATION
    # ----------------------------------------------------------------------------

    def backproject_cpu ( self, light_field, zslice = None ):
        """
        Using the given rayspread database, compute focused images at the
        given depths z_depths.

        Returns: A volume as a 3D numpy matrix that is indexed in [x,y,z]
        order.
        """

        rayspreads = self.rayspread_db.rayspreads
        z_coords = self.rayspread_db.z_coords

        nx = self.rayspread_db.nx
        ny = self.rayspread_db.ny
        ns = self.rayspread_db.ns
        nt = self.rayspread_db.nt
        nu = self.rayspread_db.nu
        nv = self.rayspread_db.nv

        # Create an empty volume
        volume = np.zeros((ny, nx, len(z_coords)), dtype = np.float32)

        # if zslice is given, compute projection for only that slice
        if zslice is not None:
            z_keys = [ k for k in rayspreads.keys() if k[0] == zslice ]
            rayspreads = extract( z_keys, rayspreads )

        # Accumulate data from the different z-slices, and add them into the lightfield.
        for spread_key in rayspreads:
            spread = rayspreads[spread_key]
            psf = spread[0]
            kernel_h_offset = spread[1]
            kernel_v_offset = spread[2]
            kernel_width = psf.shape[1]
            kernel_height = psf.shape[0]

            # Extract u,v,z
            z_idx = spread_key[0]
            u = spread_key[1]
            v = spread_key[2]

            # This performance shortcut avoids processing the rays
            # that are outside of the NA of the objective.  These rays
            # are vignetting and highly aberrated, and aren't so good
            # to use in any case.
            if (np.sqrt(np.power(u-(nu-1.0)/2.0,2) + np.power(v-(nv-1.0)/2.0,2)) >= nv/2.0):
                continue

            subaperture = light_field.subaperture(u, v)

            # zero-pad sub-aperture by half kernel width/height (otherwise splats will be cut off at
            # the left and top sizes!)
            ypad = kernel_height / 2  # setting these to 0 should disable zero padding
            xpad = kernel_width / 2 
            subaperture_pad = np.zeros((nt + 2*ypad, ns + 2*xpad))
            subaperture_pad[ ypad:ypad+nt , xpad:xpad+ns ] = subaperture

            output = filters.convolve(subaperture_pad, psf, mode='constant')

            # print "Shape: ", output.shape
            # print "Offsets: ", kernel_h_offset, kernel_v_offset
            # print "Dims: ", kernel_width, kernel_height

            # Numpy NDimage convolution crops the result
            # automatically, which removes half the kernel width and
            # height from the result.  We must adjust our kernel h and
            # v offset here to crop the output to the correct
            # dimensions.
            kernel_h_offset += kernel_width / 2
            kernel_v_offset += kernel_height / 2

            # Shift the result by the correct offset amount
            t_min = max(-kernel_v_offset,0)
            t_max = min(-kernel_v_offset+output.shape[0], output.shape[0])
            s_min = max(-kernel_h_offset,0)
            s_max = min(-kernel_h_offset+output.shape[1], output.shape[1])
            cropped_output = output[t_min:t_max, s_min:s_max]

            # Compute the opposite shift.  This is where we will sum the image in the output volume slice.
            t_min = max(kernel_v_offset,0)
            t_max = min(kernel_v_offset+output.shape[0], output.shape[0])
            s_min = max(kernel_h_offset,0)
            s_max = min(kernel_h_offset+output.shape[1], output.shape[1])
            try:
                volume_pad = np.zeros(subaperture_pad.shape)
                volume_pad[t_min:t_max, s_min:s_max] = cropped_output

                # undo zero padding
                volume[:,:,z_idx] += volume_pad[ ypad:ypad+nt , xpad:xpad+ns ]
            except ValueError:
                # Should only occur when the convolution places the
                # image out of bounds, which can happen at very high
                # NA.
                pass

        return volume

    def project_cpu( self, volume, zslice = None ):
        """
        Using the given rayspread database, compute a light field by
        projecting the volume using rayspreads.

        Option zslice allows computing the light field corresponding 
        to a single z plane in the volume.

        Returns: A lightfield as a tiled set of sub-aperture
        (i.e. pinhole) images.
        """

        rayspreads = self.rayspread_db.flipped_rayspreads
        z_coords = self.rayspread_db.z_coords

        nx = self.rayspread_db.nx
        ny = self.rayspread_db.ny
        nu = self.rayspread_db.nu
        nv = self.rayspread_db.nv

        # Create an empty volume
        lf = np.zeros((self.rayspread_db.nv * ny, self.rayspread_db.nu * nx), dtype = np.float32)

        # if zslice is given, compute projection for only that slice
        if zslice is not None:
            z_keys = [ k for k in rayspreads.keys() if k[0] == zslice ]
            rayspreads = extract( z_keys, rayspreads )

        # Accumulate data from the different z-slices, and add them into the lightfield.
        for spread_key, spread in rayspreads.items():
            # Extract u,v,z
            z_idx = spread_key[0]
            u = spread_key[1]
            v = spread_key[2]

            # This performance shortcut avoids processing the rays
            # that are outside of the NA of the objective.  These rays
            # are vignetting and highly aberrated, and aren't so good
            # to use in any case.
            if (np.sqrt(np.power(u-(nu-1.0)/2.0,2) + np.power(v-(nv-1.0)/2.0,2)) >= nv/2.0):
                continue

            kernel_h_offset = spread[1]
            kernel_v_offset = spread[2]
            psf = spread[0]
            kernel_width = psf.shape[1]
            kernel_height = psf.shape[0]
            
            zslice = volume[:,:,z_idx]

            # zero-pad zslice by half kernel width/height (otherwise
            # splats will be cut off at the left and top sizes!)
            ypad = kernel_height / 2  # setting these to 0 should disable zero padding
            xpad = kernel_width / 2 
            
            zslice_pad = np.zeros((ny + ypad, nx + xpad))
            zslice_pad[ ypad: , xpad: ] = zslice
            
            output = filters.convolve(zslice_pad, psf, mode='constant')        

            # Numpy NDimage convolution crops the result
            # automatically, which removes half the kernel width and
            # height from the result.  We must adjust our kernel h and
            # v offset here to crop the output to the correct
            # dimensions.
            kernel_h_offset += kernel_width / 2
            kernel_v_offset += kernel_height / 2

            # Shift the result by the correct offset amount
            v_min = max(-kernel_v_offset,0)
            v_max = min(-kernel_v_offset+output.shape[0], output.shape[0])
            u_min = max(-kernel_h_offset,0)
            u_max = min(-kernel_h_offset+output.shape[1], output.shape[1])
            cropped_output = output[v_min:v_max, u_min:u_max]

            # Compute the opposite shift.  This is where we will sum the image in the sub-aperture.
            v_min = max(kernel_v_offset,0)
            v_max = min(kernel_v_offset+output.shape[0], output.shape[0])
            u_min = max(kernel_h_offset,0)
            u_max = min(kernel_h_offset+output.shape[1], output.shape[1])

            
            try:
                lf_subap_pad = np.zeros(zslice_pad.shape)
                lf_subap_pad[v_min:v_max, u_min:u_max] += cropped_output

                # undo zero padding
                lf_subap = lf_subap_pad[ ypad: , xpad: ]
            
                lf[v*ny:(v+1)*ny, u*nx:(u+1)*nx] += lf_subap
            except ValueError:
                # Should only occur when the convolution places the
                # image out of bounds, which can happen at very high
                # NA.
                pass

        return LightField(lf, self.rayspread_db.nu, self.rayspread_db.nv,
                          self.rayspread_db.ns, self.rayspread_db.nt,
                          representation = LightField.TILED_SUBAPERTURE)


    # ----------------------------------------------------------------------------
    #                          SART METHODS (DEPRECATED)
    # ----------------------------------------------------------------------------
    
    def backproject_angle ( self, light_field_angle, u, v):
        """
        Using the given rayspread database and SINGLE lightfield pinhole
        image (for SART, not SIRT), backproject into volume.

        Returns: A volume as a 3D numpy matrix that is indexed in [x,y,z]
        order.
        """
        rayspreads = self.rayspread_db.rayspreads
        z_coords = self.rayspread_db.z_coords

        nx = self.rayspread_db.ns
        ny = self.rayspread_db.nt

        # Create an empty volume
        volume = np.zeros((ny, nx, len(z_coords)), dtype = np.float32)

        # Accumulate data from the different z-slices, and add them into the lightfield.
        for z in range(len(z_coords)):  
            spread_key=(z,u,v)
            if spread_key in rayspreads:
                spread = rayspreads[spread_key]
                psf = spread[0]
                kernel_h_offset = spread[1]
                kernel_v_offset = spread[2]
                kernel_width = psf.shape[1]
                kernel_height = psf.shape[0]

                # Extract u,v,z
                z_idx = spread_key[0]
                u = spread_key[1]
                v = spread_key[2]

                subaperture = light_field_angle
                output = filters.convolve(subaperture, psf, mode='constant')

                # Numpy NDimage convolution shifts the image by half the kernel width,
                # so we must compensate for that offset here.
                kernel_h_offset += kernel_width / 2
                kernel_v_offset += kernel_height / 2

                # Shift the result by the correct offset amount
                t_min = max(-kernel_v_offset,0)
                t_max = min(-kernel_v_offset+output.shape[0], output.shape[0])
                s_min = max(-kernel_h_offset,0)
                s_max = min(-kernel_h_offset+output.shape[1], output.shape[1])
                cropped_output = output[t_min:t_max, s_min:s_max]

                # Compute the opposite shift.  This is where we will sum the image in the output volume slice.
                t_min = max(kernel_v_offset,0)
                t_max = min(kernel_v_offset+output.shape[0], output.shape[0])
                s_min = max(kernel_h_offset,0)
                s_max = min(kernel_h_offset+output.shape[1], output.shape[1])
                volume[t_min:t_max, s_min:s_max, z_idx] += cropped_output

        return volume


    def project_angle( self, volume, u, v ):
        """
        Using the given rayspread database, compute a SINGLE light field view by
        projecting the volume using rayspreads (for SART, not SIRT).

        Returns: A single lightfield view (pinhole image)
        """

        rayspreads = self.rayspread_db.flipped_rayspreads
        z_coords = self.rayspread_db.z_coords

        nx = volume.shape[1]
        ny = volume.shape[0]

        # Create an empty lightfield pinhole image
        lf = np.zeros((ny, nx), dtype = np.float32)

        # Accumulate data from the different z-slices, and add them into the lightfield pinhole image.
        for z in range(len(z_coords)):  
            spread_key=(z,u,v)
            if spread_key in rayspreads:
                spread = rayspreads[spread_key]
                psf = spread[0]
                kernel_h_offset = spread[1]
                kernel_v_offset = spread[2]
                kernel_width = psf.shape[1]
                kernel_height = psf.shape[0]

                zslice = volume[:,:,z]
                output = filters.convolve(zslice, psf, mode='constant')

                # Numpy NDimage convolution shifts the image by half the kernel width,
                # so we must compensate for that offset here.
                kernel_h_offset += kernel_width / 2
                kernel_v_offset += kernel_height / 2

                # Shift the result by the correct offset amount
                v_min = max(-kernel_v_offset,0)
                v_max = min(-kernel_v_offset+output.shape[0], output.shape[0])
                u_min = max(-kernel_h_offset,0)
                u_max = min(-kernel_h_offset+output.shape[1], output.shape[1])
                cropped_output = output[v_min:v_max, u_min:u_max]

                # Compute the opposite shift.  This is where we will sum the image in the sub-aperture.
                v_min = max(kernel_v_offset,0)
                v_max = min(kernel_v_offset+output.shape[0], output.shape[0])
                u_min = max(kernel_h_offset,0)
                u_max = min(kernel_h_offset+output.shape[1], output.shape[1])

                lf[v_min:v_max, u_min:u_max] += cropped_output

        return lf 

    # -----------------------------------------------------------------------------
    #                      MISC. FUNCTIONS FOR DEBUGGING
    # -----------------------------------------------------------------------------

    def compare_cpu_gpu_performance (self):
        '''
        Compare the performance of the GPU and CPU implementations of
        the forward and back projection algorithms.
        '''

        print 'Testing GPU vs. CPU SIRT performance...'

        # Create a synthetic light field and volume for testing purposes
        im = np.ones((self.rayspread_db.nt*self.rayspread_db.nv,
                      self.rayspread_db.ns*self.rayspread_db.nu), dtype=np.float32)
        lf = LightField(im, self.rayspread_db.nu, self.rayspread_db.nv,
                        self.rayspread_db.ns, self.rayspread_db.nt,
                        representation = LightField.TILED_SUBAPERTURE)

        vol = np.ones((self.rayspread_db.ny,
                       self.rayspread_db.nx,
                       self.rayspread_db.nz), dtype=np.float32)

        print '\t--> Forward projection:'
        import time
        t0 = time.time()
        self.project_gpu(vol)
        gpu_dt = time.time()-t0
        t0 = time.time()
        self.project_cpu(vol)
        cpu_dt = time.time()-t0
        print '\t    GPU: %0.2fs  CPU: %0.2fs   Speedup: %0.1fx' % (gpu_dt, cpu_dt, cpu_dt/gpu_dt)

        print '\t--> Back projection:'
        import time
        t0 = time.time()
        vol = self.backproject_gpu(lf)
        gpu_dt = time.time()-t0

        t0 = time.time()
        self.backproject_cpu(lf)
        cpu_dt = time.time()-t0
        print '\t    GPU: %0.2fs  CPU: %0.2fs   Speedup: %0.1fx' % (gpu_dt, cpu_dt, cpu_dt/gpu_dt)



    def test_forward_project(self):
        '''
        Compares the acclerated GPU implemenation to the CPU implementation to
        verify that they return the same results, modulo floating point error.
        '''
        print 'testing SIRT project'
        vol_zero = np.zeros((self.rayspread_db.ny,
                             self.rayspread_db.nx,
                             self.rayspread_db.nz), dtype=np.float32)
        for z in range(self.rayspread_db.nz):
            # generate volume with a single '1' voxel at position (x,y,z)
            delta_func = np.copy(vol_zero)
            x = int(self.rayspread_db.nx/2)
            y = int(self.rayspread_db.ny/2)
            delta_func[y, x, z] = 1
            
            # forward project with CPU code only!
            lf_cpu = self.project_cpu(delta_func)
            lf_gpu = self.project_gpu(delta_func)

            lf_gpu_im = lf_gpu.asimage(representation = LightField.TILED_LENSLET)
            lf_cpu_im = lf_cpu.asimage(representation = LightField.TILED_LENSLET)

            filename = 'psf/lf_cpu_%d.tif' % (z)
            save_image(filename, lf_cpu_im, dtype=np.float32)
            filename = 'psf/lf_gpu_%d.tif' % (z)
            save_image(filename, lf_gpu_im, dtype=np.float32)
            
            # compare the results to see if they match
            error = np.abs(lf_cpu_im - lf_gpu_im)
            print 'lf totalerr:', np.sum(error), 'maxerr:', error.max(), 'cpu/gpu .max():', lf_cpu_im.max(), lf_gpu_im.max()

    def test_back_project(self):
        '''
        Compares the acclerated GPU implemenation to the CPU implementation to
        verify that they return the same results, modulo floating point error.
        '''

        # check that CPU and GPU backproject() functions are working
        # (this code doesn't check project(), only backproject() )            
        print 'testing SIRT back project'
        vol_zero = np.zeros((self.rayspread_db.ny,
                             self.rayspread_db.nx,
                             self.rayspread_db.nz), dtype=np.float32)
            
        for z in range(self.rayspread_db.nz):
            # generate volume with a single '1' voxel at position (x,y,z)
            delta_func = np.copy(vol_zero)
            x = int(self.rayspread_db.nx/2)
            y = int(self.rayspread_db.ny/2)
            delta_func[y, x, z] = 1

            # forward project with CPU code only!
            lf_cpu = self.project_cpu(delta_func)

            # back project with CPU and GPU code                
            v_cpu = self.backproject_cpu(lf_cpu)
            v_gpu = self.backproject_gpu(lf_cpu)

            # compare the results to see if they match
            error = np.abs(v_cpu - v_gpu)
            print 'v total error:', np.sum(error), 'max error:', error.max(), 'cpu/gpu .max():',    v_cpu.max(), v_gpu.max()
                
    def test_forward_backward(self):
        '''
        Tests to make sure that a voxel that is forward projected then
        backprojected ends up in the same place.
        '''
        vol_zero = self.rayspread_db.empty_volume()
        
        count = 0 # keep track of how many z-planes are correct
        for z in range(self.rayspread_db.nz):
#        for z in range(9,10):
            for x_offset in range(self.rayspread_db.supersample_factor):
                for y_offset in range(self.rayspread_db.supersample_factor):
                    print '\t--> z = ', z, ' x_o = ', x_offset, ' y_o = ', y_offset

                    # generate volume with a single '1' voxel at position (x,y,z)
                    delta_func = np.copy(vol_zero)
                    x = int(self.rayspread_db.nx/2 + x_offset)
                    y = int(self.rayspread_db.ny/2 + y_offset)
                    delta_func[y, x, z] = 1

                    # forward and back project, to get a focal stack
                    lf_gpu = self.project_gpu(delta_func)

                    filename = 'psf/lf_%d_%d_%d.tif' % (z, x_offset, y_offset)
                    save_image(filename,
                               lf_gpu.asimage(representation = LightField.TILED_LENSLET), dtype=np.float32)
                    
                    v_gpu = self.backproject_gpu(lf_gpu)

                    filename = 'psf/vol_%d_%d_%d.tif' % (z, x_offset, y_offset)
                    save_image(filename, v_gpu, dtype=np.float32)

                    # get indices (y',x',z') of max voxel of focal stack
                    v_index =  np.unravel_index(np.argmax(v_gpu), v_gpu.shape)
                    delta_index =  np.unravel_index(np.argmax(delta_func), v_gpu.shape)

                    # check if indices (x,y,z) match (x',y',z')--they should if everything's right!
                    if np.argmax(v_gpu) - np.argmax(delta_func) == 0:
                        print '\t    index of maximum matches', delta_index
                        count += 1
                    else:
                        print '\t    index of maximum does not match!', 'focal stack:', v_index,'original:', delta_index

        print 'position of max intensity voxel matched input for', count, 'of', self.rayspread_db.nz, 'z-planes'
        print 'finished testing SIRT projection functions, quitting...'

