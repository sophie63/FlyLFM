# __BEGIN_LICENSE__
#
# Copyright 2012 Stanford University.  All rights reserved.
#
# __END_LICENSE__

'''
This file contains various utilities and routines that are useful for
direct deconvolution methods. Note that some of these routines are
only useful for debugging or solving small toy problems.
'''
import sys, time
import numpy as np
from lflib.volume import LightFieldProjection 
from lflib.imageio import save_image
from lflib.volume import roundUp, padVector
from scipy.sparse import lil_matrix, coo_matrix, eye, vstack, hstack, identity, triu

try:
    import pyopencl as cl
    import pyopencl.array as cl_array
    import pyopencl.characterize as cl_characterize
    LFLIB_HAVE_OPENCL = True
except ImportError:
    LFLIB_HAVE_OPENCL = False


# ------------------------------------------------------------------------------
#                               OPENCL KERNELS
# ------------------------------------------------------------------------------

src = '''

__kernel void compute_submatrix_reference(__global float* output,
                                __global int* splat1_x_coords, __global int* splat1_y_coords, 
                                __global int* splat2_x_coords, __global int* splat2_y_coords,
                                __global float* splat1_vals, __global float* splat2_vals,
                                int n_coords1, int n_coords2, int nu, int nv, int nx, int ny, int yy,
                                __local int* scratch_x, __local int* scratch_y, 
                                __local float* scratch_vals) {

  // Copy frequently used data to local memory
  //  for (int j = 0; j < n_coords2/4; ++j) {
  //    ((__local int4*)scratch_x)[j] = ((__global int4*)splat2_x_coords)[j];
  //    ((__local int4*)scratch_y)[j] = ((__global int4*)splat2_y_coords)[j];
  //    ((__local float4*)scratch_vals)[j] = ((__global float4*)splat2_vals)[j];
  //  }
  //  barrier(CLK_LOCAL_MEM_FENCE);

  // Only do data processing for valid voxels
  if ((int)get_global_id(0) >= ny || (int)get_global_id(1) >= nx)
    return;

  // Which row of the covariance block are we working on?  
  int row_idx = get_global_id(0)*nx + get_global_id(1);

  // Each work unit operates on one row of the covariance block.  The
  // ... offset varies depending on which row we are working on.
  // (note that centering could be off, causing asymmetric edge effects)
  int ry = get_global_id(0) - ny/2; int rx = get_global_id(1) - nx/2;

  // We iterate over the columns of the submatrix block here.  The 
  // offset of splat2 is determined based on which column we are working
  // on.
  //  
  int dy = yy - ny/2;
  for (int xx = 0; xx < nx; ++xx) {
    int dx = xx - nx/2;
    int col_idx = xx;

    // skip redundant lower-triangular entries.  We will fill those in later.
    if (yy*nx + col_idx < row_idx)
      continue;

    // Finally we compute the inner product for this particular
    // offset of splat1 and splat2.  The inner product can be computed
    // efficiently by a small number of comparisons between of the
    // coordinates of nonzero entries in splat1 and spalt2.  This is only
    // fast if the splats are very sparse.
    //
    float sum = 0;
    for (int i = 0; i < n_coords1; ++i) {
      int coords1_x = splat1_x_coords[i]+(rx*nu);
      int coords1_y = splat1_y_coords[i]+(ry*nv);

      int in_bounds1 = (coords1_x >= 0) && (coords1_y >= 0) && (coords1_x < nx*nu) && (coords1_y < ny*nv);
      float val1 = splat1_vals[i];
      
      for (int j = 0; j < n_coords2; ++j) {
         int coords2_x = splat2_x_coords[j]+(dx*nu);
         int coords2_y = splat2_y_coords[j]+(dy*nv);

         int coords_match = (coords2_x == coords1_x) && (coords2_y == coords1_y);
         int in_bounds2 = (coords2_x >= 0) && (coords2_y >= 0) && (coords2_x < nx*nu) && (coords2_y < ny*nv);
         float val2 = splat2_vals[j];
        
         // These are all boolean tests to see if the coords overlap and if they are
         // in the bounds that we will consider for the sum.  This could be an if statement,
         // but we will avoid branching in the GPU kernel and make this a series of comparisons
         // and multiplications instead.  The effect is the same.
         sum += in_bounds1 * in_bounds2 * coords_match * val1 * val2;
      }
    }
    output[row_idx*nx+col_idx] = sum;
  }
}

__kernel void compute_submatrix(__global float* output,
                                __global int* splat1_x_coords, __global int* splat1_y_coords, 
                                __global int* splat2_x_coords, __global int* splat2_y_coords,
                                __global float* splat1_vals, __global float* splat2_vals,
                                int n_coords1, int n_coords2, int nu, int nv, int nx, int ny, int yy ) { //,
                                //__local int* scratch_x,
                                //__local int* scratch_y, 
                                //__local float* scratch_vals) {

  // Copy frequently used data to local memory
  //for (int j = 0; j < n_coords2/4; ++j) {
  //  ((__local int4*)scratch_x)[j] = ((__global int4*)splat2_x_coords)[j];
  //  ((__local int4*)scratch_y)[j] = ((__global int4*)splat2_y_coords)[j];
  //  ((__local float4*)scratch_vals)[j] = ((__global float4*)splat2_vals)[j];
  //}

  //barrier(CLK_LOCAL_MEM_FENCE);

  // Only do data processing for valid voxels
  if ((int)get_global_id(0) >= ny || (int)get_global_id(1) >= nx)
    return;

  // Which row of the covariance block are we working on?  
  int row_idx = get_global_id(0)*nx + get_global_id(1);

  // Each work unit operates on one row of the covariance block.  The
  // offset varies depending on which row we are working on.
  // (note that centering here could be slightly off, causing asymmetric edge effects)
  int ry = get_global_id(0) - ny/2; int rx = get_global_id(1) - nx/2;

  // We iterate over the columns of the submatrix block here.  The 
  // offset of splat2 is determined based on which column we are working
  // on.
  // yy is the index of the block of columns of the covariance
  // matrix block that we are currently computing. 
  // xx iterates over the entries in a particular column.
  //
  int dy = yy - ny/2; 
  for (int xx = 0; xx < nx; ++xx) {
    int dx = xx - nx/2;
    int col_idx = xx;

    // skip redundant lower-triangular entries 
    if (yy*nx + col_idx < row_idx)
      continue;

    // Finally we compute the inner product for this particular
    // offset of splat1 and splat2.  The inner product can be computed
    // efficiently by a small number of comparisons between of the
    // coordinates of nonzero entries in splat1 and spalt2.  This is only
    // fast if the splats are very sparse.
    //
    float sum = 0.0;
    for (int i = 0; i < n_coords1/4; ++i) {
      
      int4 coords1_x = ((__global int4*)splat1_x_coords)[i]+(rx*nu);
      int4 coords1_y = ((__global int4*)splat1_y_coords)[i]+(ry*nv);

      int4 in_bounds1 = ((coords1_x >= 0) && (coords1_y >= 0) && (coords1_x < nx*nu) && (coords1_y < ny*nv));
      float4 val1 = ((__global float4*)splat1_vals)[i];

      for (int j = 0; j < n_coords2/4; ++j) {
        int4 coords2_x = ((__global int4*)splat2_x_coords)[j]+(dx*nu);
        int4 coords2_y = ((__global int4*)splat2_y_coords)[j]+(dy*nv);

        int4 in_bounds2 = ((coords2_x >= 0) && (coords2_y >= 0) && (coords2_x < nx*nu) && (coords2_y < ny*nv));
        float4 val2 = ((__global float4*)splat2_vals)[j];
        
        // These are all boolean tests to see if the coords overlap and if they are
        // in the bounds that we will consider for the sum.  This could be an if statement,
        // but we will avoid branching in the GPU kernel and make this a series of comparisons
        // and multiplications instead.  The effect is the same.

        sum += -in_bounds1.x * -in_bounds2.x * -(coords1_x.x == coords2_x.x) * -(coords1_y.x == coords2_y.x) * val1.x * val2.x;
        sum += -in_bounds1.x * -in_bounds2.y * -(coords1_x.x == coords2_x.y) * -(coords1_y.x == coords2_y.y) * val1.x * val2.y;
        sum += -in_bounds1.x * -in_bounds2.z * -(coords1_x.x == coords2_x.z) * -(coords1_y.x == coords2_y.z) * val1.x * val2.z;
        sum += -in_bounds1.x * -in_bounds2.w * -(coords1_x.x == coords2_x.w) * -(coords1_y.x == coords2_y.w) * val1.x * val2.w;

        sum += -in_bounds1.y * -in_bounds2.x * -(coords1_x.y == coords2_x.x) * -(coords1_y.y == coords2_y.x) * val1.y * val2.x;
        sum += -in_bounds1.y * -in_bounds2.y * -(coords1_x.y == coords2_x.y) * -(coords1_y.y == coords2_y.y) * val1.y * val2.y;
        sum += -in_bounds1.y * -in_bounds2.z * -(coords1_x.y == coords2_x.z) * -(coords1_y.y == coords2_y.z) * val1.y * val2.z;
        sum += -in_bounds1.y * -in_bounds2.w * -(coords1_x.y == coords2_x.w) * -(coords1_y.y == coords2_y.w) * val1.y * val2.w;

        sum += -in_bounds1.z * -in_bounds2.x * -(coords1_x.z == coords2_x.x) * -(coords1_y.z == coords2_y.x) * val1.z * val2.x;
        sum += -in_bounds1.z * -in_bounds2.y * -(coords1_x.z == coords2_x.y) * -(coords1_y.z == coords2_y.y) * val1.z * val2.y;
        sum += -in_bounds1.z * -in_bounds2.z * -(coords1_x.z == coords2_x.z) * -(coords1_y.z == coords2_y.z) * val1.z * val2.z;
        sum += -in_bounds1.z * -in_bounds2.w * -(coords1_x.z == coords2_x.w) * -(coords1_y.z == coords2_y.w) * val1.z * val2.w;

        sum += -in_bounds1.w * -in_bounds2.x * -(coords1_x.w == coords2_x.x) * -(coords1_y.w == coords2_y.x) * val1.w * val2.x;
        sum += -in_bounds1.w * -in_bounds2.y * -(coords1_x.w == coords2_x.y) * -(coords1_y.w == coords2_y.y) * val1.w * val2.y;
        sum += -in_bounds1.w * -in_bounds2.z * -(coords1_x.w == coords2_x.z) * -(coords1_y.w == coords2_y.z) * val1.w * val2.z;
        sum += -in_bounds1.w * -in_bounds2.w * -(coords1_x.w == coords2_x.w) * -(coords1_y.w == coords2_y.w) * val1.w * val2.w;
      }
    }
    output[row_idx*nx+col_idx] = sum;
  }
}

__kernel void compute_submatrix_fast(__global float* output,
                                __global int* splat1_x_coords, __global int* splat1_y_coords, 
                                __global int* splat2_x_coords, __global int* splat2_y_coords,
                                __global float* splat1_vals, __global float* splat2_vals,
                                int n_coords1, int n_coords2, int nu, int nv, int nx, int ny, int yy ) { //,
                                //__local int* scratch_x,
                                //__local int* scratch_y, 
                                //__local float* scratch_vals) {

  
  // Only do data processing for valid voxels
  if ((int)get_global_id(0) >= ny || (int)get_global_id(1) >= nx)
    return;

  // Which row of the covariance block are we working on?  
  int row_idx = get_global_id(0)*nx + get_global_id(1);

  // Each work unit operates on one row of the covariance block.  The
  // offset varies depending on which row we are working on.
  // (note that centering here could be slightly off, causing asymmetric edge effects)
  int ry = get_global_id(0) - ny/2; int rx = get_global_id(1) - nx/2;

  // We iterate over the columns of the submatrix block here.  The 
  // offset of splat2 is determined based on which column we are working
  // on.
  // yy is the index of the block of columns of the covariance
  // matrix block that we are currently computing. 
  // xx iterates over the entries in a particular column.
  //
  int dy = yy - ny/2; 
  for (int xx = 0; xx < nx; ++xx) {
    int dx = xx - nx/2;
    int col_idx = xx;

    // skip redundant lower-triangular entries 
    //if (yy*nx + col_idx < row_idx)
    //  continue;

    // Finally we compute the inner product for this particular
    // offset of splat1 and splat2.  The inner product can be computed
    // efficiently by a small number of comparisons between of the
    // coordinates of nonzero entries in splat1 and spalt2.  This is only
    // fast if the splats are very sparse.
    //
    float sum = 0.0;
    for (int i = 0; i < n_coords1/4; ++i) {
      
      int4 coords1_x = ((__global int4*)splat1_x_coords)[i]+(rx*nu);
      int4 coords1_y = ((__global int4*)splat1_y_coords)[i]+(ry*nv);

      int4 in_bounds1 = (1 == 1); //((coords1_x >= 0) && (coords1_y >= 0) && (coords1_x < nx*nu) && (coords1_y < ny*nv));
      float4 val1 = ((__global float4*)splat1_vals)[i];

      for (int j = 0; j < n_coords2/4; ++j) {
        int4 coords2_x = ((__global int4*)splat2_x_coords)[j]+(dx*nu);
        int4 coords2_y = ((__global int4*)splat2_y_coords)[j]+(dy*nv);

        int4 in_bounds2 = (1 == 1); //((coords2_x >= 0) && (coords2_y >= 0) && (coords2_x < nx*nu) && (coords2_y < ny*nv));
        float4 val2 = ((__global float4*)splat2_vals)[j];
        
        // These are all boolean tests to see if the coords overlap and if they are
        // in the bounds that we will consider for the sum.  This could be an if statement,
        // but we will avoid branching in the GPU kernel and make this a series of comparisons
        // and multiplications instead.  The effect is the same.

        sum += -in_bounds1.x * -in_bounds2.x * -(coords1_x.x == coords2_x.x) * -(coords1_y.x == coords2_y.x) * val1.x * val2.x;
        sum += -in_bounds1.x * -in_bounds2.y * -(coords1_x.x == coords2_x.y) * -(coords1_y.x == coords2_y.y) * val1.x * val2.y;
        sum += -in_bounds1.x * -in_bounds2.z * -(coords1_x.x == coords2_x.z) * -(coords1_y.x == coords2_y.z) * val1.x * val2.z;
        sum += -in_bounds1.x * -in_bounds2.w * -(coords1_x.x == coords2_x.w) * -(coords1_y.x == coords2_y.w) * val1.x * val2.w;

        sum += -in_bounds1.y * -in_bounds2.x * -(coords1_x.y == coords2_x.x) * -(coords1_y.y == coords2_y.x) * val1.y * val2.x;
        sum += -in_bounds1.y * -in_bounds2.y * -(coords1_x.y == coords2_x.y) * -(coords1_y.y == coords2_y.y) * val1.y * val2.y;
        sum += -in_bounds1.y * -in_bounds2.z * -(coords1_x.y == coords2_x.z) * -(coords1_y.y == coords2_y.z) * val1.y * val2.z;
        sum += -in_bounds1.y * -in_bounds2.w * -(coords1_x.y == coords2_x.w) * -(coords1_y.y == coords2_y.w) * val1.y * val2.w;

        sum += -in_bounds1.z * -in_bounds2.x * -(coords1_x.z == coords2_x.x) * -(coords1_y.z == coords2_y.x) * val1.z * val2.x;
        sum += -in_bounds1.z * -in_bounds2.y * -(coords1_x.z == coords2_x.y) * -(coords1_y.z == coords2_y.y) * val1.z * val2.y;
        sum += -in_bounds1.z * -in_bounds2.z * -(coords1_x.z == coords2_x.z) * -(coords1_y.z == coords2_y.z) * val1.z * val2.z;
        sum += -in_bounds1.z * -in_bounds2.w * -(coords1_x.z == coords2_x.w) * -(coords1_y.z == coords2_y.w) * val1.z * val2.w;

        sum += -in_bounds1.w * -in_bounds2.x * -(coords1_x.w == coords2_x.x) * -(coords1_y.w == coords2_y.x) * val1.w * val2.x;
        sum += -in_bounds1.w * -in_bounds2.y * -(coords1_x.w == coords2_x.y) * -(coords1_y.w == coords2_y.y) * val1.w * val2.y;
        sum += -in_bounds1.w * -in_bounds2.z * -(coords1_x.w == coords2_x.z) * -(coords1_y.w == coords2_y.z) * val1.w * val2.z;
        sum += -in_bounds1.w * -in_bounds2.w * -(coords1_x.w == coords2_x.w) * -(coords1_y.w == coords2_y.w) * val1.w * val2.w;
      }
    }
    output[row_idx*nx+col_idx] = sum;
  }
}

'''

# ------------------------------------------------------------------------------
#                               UTILITY / DEBUGGING
# ------------------------------------------------------------------------------
    
class DirectReconstruction(LightFieldProjection):

    def __init__(self, rayspread_db, disable_gpu = False, gpu_id = None):

        self.rayspread_db = rayspread_db
        self.disable_gpu = disable_gpu
        self.gpu_id = gpu_id

        # Initialize base class
        LightFieldProjection.__init__(self, rayspread_db)
        
        if LFLIB_HAVE_OPENCL:
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

            self.cl_prog = cl.Program(self.cl_ctx, src).build()
            self.cl_queue = cl.CommandQueue(self.cl_ctx)

            # Set up OpenCL
            self.WORKGROUP_X = 16
            self.WORKGROUP_Y = 16
            self.WORKGROUP_Z = 1

        else:
            print "ERROR: direct decovolution require OpenCL."
            raise SystemExit

    # COVARIANCE MATRIX
    #
    # These routines compute the covariance matrix transpose(A)*A.
    # This matrix is much smaller than A itself (the covariance matrix
    # has Nvoxel x Nvoxel entires), but it is nonetheless useful for
    # solving using ADMM and other direct algorithms.
    #
    def covA_submatrix( self, sirt, z1, z2, reference_kern = False, fast_method = True ):
        nx = self.rayspread_db.nx
        ny = self.rayspread_db.ny
        nu = self.rayspread_db.nu
        nv = self.rayspread_db.nv

        lf1 = sirt.splat( (self.rayspread_db.ny/2,self.rayspread_db.nx/2,z1) ).asimage().astype(np.float32)
        lf2 = sirt.splat( (self.rayspread_db.ny/2,self.rayspread_db.nx/2,z2) ).asimage().astype(np.float32)

        nonzero1 = np.nonzero(lf1)
        nonzero2 = np.nonzero(lf2)
        splat1_y_coords_buf = cl_array.to_device(self.cl_queue, padVector(nonzero1[0],4).astype(np.int32))
        splat1_x_coords_buf = cl_array.to_device(self.cl_queue, padVector(nonzero1[1],4).astype(np.int32))
        splat2_y_coords_buf = cl_array.to_device(self.cl_queue, padVector(nonzero2[0],4).astype(np.int32))
        splat2_x_coords_buf = cl_array.to_device(self.cl_queue, padVector(nonzero2[1],4).astype(np.int32))

        vals1 = padVector(lf1[nonzero1],4)
        vals2 = padVector(lf2[nonzero2],4)
        n_coords1 = len(vals1)
        n_coords2 = len(vals2)
        splat1_vals_buf = cl_array.to_device(self.cl_queue, vals1)
        splat2_vals_buf = cl_array.to_device(self.cl_queue, vals2)

        # get radius of splat overlap
        #radius = (np.sqrt( (np.max(nonzero1[0]) - lf_center[0])**2 + (np.max(nonzero1[1]) - lf_center[1])**2 ) + 
        #np.sqrt( (np.max(nonzero2[0]) - lf_center[0])**2 + (np.max(nonzero2[1]) - lf_center[1])**2 ) )

        localSize = (self.WORKGROUP_X, self.WORKGROUP_Y)
        globalSize = (roundUp(ny,self.WORKGROUP_Y),
                      roundUp(nx,self.WORKGROUP_X))

        if reference_kern:
            kern = self.cl_prog.compute_submatrix_reference
        else:
            if fast_method:
                kern = self.cl_prog.compute_submatrix_fast
            else:
                kern = self.cl_prog.compute_submatrix

        submatrix_buf = cl_array.zeros(self.cl_queue, (nx*ny,nx), dtype=np.float32)
        kern.set_arg(0, submatrix_buf.data)
        
        kern.set_arg(1, splat1_x_coords_buf.data)
        kern.set_arg(2, splat1_y_coords_buf.data)
        kern.set_arg(3, splat2_x_coords_buf.data)
        kern.set_arg(4, splat2_y_coords_buf.data)
        kern.set_arg(5, splat1_vals_buf.data)
        kern.set_arg(6, splat2_vals_buf.data)
        kern.set_arg(7, np.int32(n_coords1))
        kern.set_arg(8, np.int32(n_coords2))
        kern.set_arg(9, np.int32(nu))
        kern.set_arg(10, np.int32(nv))
        kern.set_arg(11, np.int32(nx))
        kern.set_arg(12, np.int32(ny))
        #kern.set_arg(14, cl.LocalMemory(n_coords2*4))
        #kern.set_arg(15, cl.LocalMemory(n_coords2*4))
        #kern.set_arg(16, cl.LocalMemory(n_coords2*4))

        # Execute the kernel
        # yy will be the index passed to the gpu kernel 
        # indicating which set of nx columns will be computed
        submatrix = lil_matrix( (nx*ny, nx*ny) )
        if not fast_method:
            for yy in range(ny):
                # print yy, nx, ny
                submatrix_buf.fill(0.0, self.cl_queue)
                kern.set_arg(13, np.int32(yy))
                #t0 = time.time()
                cl.enqueue_nd_range_kernel(self.cl_queue, kern, globalSize, localSize)
                submatrix[:,yy*nx:yy*nx+nx] = submatrix_buf.get()
                # print 'elapsed time = ', time.time() - t0
        else:
            print ' using fast compute_submatrix() method, which ignores edge effects'
            
            # first compute first nx columns (we actually only need 1 column)
            # this can be optimized...there may also be a bug with supersampling at native plane
            yy = 0
            submatrix_buf.fill(0.0, self.cl_queue)
            kern.set_arg(13, np.int32(yy))
            #t0 = time.time()
            cl.enqueue_nd_range_kernel(self.cl_queue, kern, globalSize, localSize)
            submatrix[:,yy*nx:yy*nx+nx] =submatrix_buf.get()

            print 'now filling in the rest of the matrix'
            # given one column, fill in the rest of the columns
            # this can be optimized...
            for y1 in range(ny):
                for x1 in range(nx):
                    for y2 in range(ny):
                        for x2 in range(nx):
                            dy = abs(y1 - y2)
                            dx = abs(x1 - x2)
                            submatrix[y2*nx + x2, y1*nx + x1] = submatrix[dy*nx + dx,0]
                
            
        self.cl_queue.finish()

        # Reflect across the main diagonal to fill in the rest of the
        # entries in the matrix.
        submatrix = triu(submatrix,1).T + triu(submatrix,0)

        # Return the result
        return submatrix
        
    def covA(self, raydb, test = False, reference_kern = False):

        import time
        t0 = time.time()
        sirt = SimultaneousIterativeReconstruction(raydb, self.disable_gpu, self.gpu_id)

        print ('Constructing covariance matrix for A with %d rows and %d cols' % 
                  (raydb.nx*raydb.ny*raydb.nz, raydb.nx*raydb.ny*raydb.nz))
        AtA = np.zeros((raydb.nx*raydb.ny*raydb.nz, raydb.nx*raydb.ny*raydb.nz));
        for z1 in range(raydb.nz):
            for z2 in range(z1, raydb.nz):
                print '\t--> Submatrix: ', z1, z2
                if reference_kern:
                    Asub = self.covA_submatrix(sirt, z1, z2, reference_kern=True)
                    filename = 'covA/covA_submatrix_reference_%d_%d.png' % (z1, z2)
                else:
                    Asub = self.covA_submatrix(sirt, z1, z2)
                    filename = 'covA/covA_submatrix_%d_%d.png' % (z1, z2)

                if test:
                    Adense = Asub.todense()
                    print '\t    ', Adense.min(), Adense.max(), Adense.mean()
                    save_image(filename, Adense/Adense.max()*65535, dtype=np.uint16)

                rowstart = z1*raydb.nx*raydb.ny
                rowend = z1*raydb.nx*raydb.ny+raydb.nx*raydb.ny
                colstart = z2*raydb.nx*raydb.ny
                colend = z2*raydb.nx*raydb.ny+raydb.nx*raydb.ny
                AtA[rowstart:rowend,colstart:colend] = Adense

        print 'Total time elapsed:', time.time()-t0,' seconds.'
        if test:
            save_image("covA/AtA_gpu.png", AtA/AtA.max()*65535, dtype=np.uint16)      

        return AtA

    # DEBUGGING ROUTINES
    #
    # These routines compute the full projection matrix ("A") or columns thereof.
    # These are rarely useful in practice, since A can be quite large.
    # However, these routines can be useful for debugging small toy
    # problems.
    #
    def projection_matrix_column(self, sirt, raydb, column):
        '''
        Compute a column of the projection matrix.  NOTE: this function is
        *only* useful for small toy problems!  The projection matrix is usually
        too large to compute!
        '''
        idx = column
        z = int(idx/(raydb.nx*raydb.ny))
        idx -= z*raydb.nx*raydb.ny
        y = int(idx/raydb.nx) 
        idx -= y*raydb.nx
        x = idx
        lf = sirt.splat( (y,x,z) )

        # For debugging: save out splats
        #        fname = 'covA/splat_%d_%d_%d.png' % (y, x, z)
        #        test = lf.asimage()
        #        save_image(fname, test/test.max()*65535, dtype=np.uint16)
        return lf.asimage().reshape(raydb.ns*raydb.nt*raydb.nu*raydb.nv)

    def projection_matrix(self, raydb, cov=False):
        '''
        Constructs the full projection matrix.  NOTE: this function is *only*
        useful for small toy problems!  The projection matrix is usually too
        large to compute or to store!
        '''
        sirt = LightFieldProjection(raydb, self.disable_gpu, self.gpu_id)
        A = np.zeros((raydb.ns*raydb.nt*raydb.nu*raydb.nv, raydb.nx*raydb.ny*raydb.nz))
        for z in range(raydb.nz):
            for y in range(raydb.ny):
                for x in range(raydb.nx):
                    column_idx = z*raydb.nx*raydb.ny+y*raydb.nx+x
                    A[:,column_idx] = self.projection_matrix_column(sirt, raydb, column_idx)

        # return projection matrix and covariance if cov is true, otherwise just the projection matrix
        if cov:
            AtA = np.dot(A.T,A)
            save_image("covA/AtA_simpledirect.png", AtA / AtA.max() * 65535, dtype=np.uint16)
            AtA = AtA[0:raydb.nx*raydb.ny,0:raydb.nx*raydb.ny]
            save_image("covA/direct_submatrix_0_0.png", AtA / AtA.max() * 65535, dtype=np.uint16)
            return A, AtA
        else:
            return A

def flatten_volume(vol):
    # flatten volume, using code from direct_deconvolution.projection_matrix_column
    # uses 1 0 2 ordering (neither row nor column major...to be changed)
    nx = vol.shape[1]
    ny = vol.shape[0]
    nz = vol.shape[2]
    vol_flat = np.zeros(nx*ny*nz)
    for idx  in range(vol_flat.shape[0]):
        id = idx
        z = int(idx/(nx*ny))
        idx -= z*nx*ny
        y = int(idx/nx)
        idx -= y*nx
        x = idx
        vol_flat[id] = vol[y,x,z]
    return vol_flat
    
def unflatten_volume(vol_flat, shape):
    # unflatten volume, using code from direct_deconvolution.projection_matrix_column
    # uses 1 0 2 ordering (neither row nor column major...to be changed)
    ny = shape[0]
    nx = shape[1]
    nz = shape[2]
    vol = np.zeros((ny,nx,nz))
    for idx  in range(vol_flat.shape[0]):
        id = idx
        z = int(idx/(nx*ny))
        idx -= z*nx*ny
        y = int(idx/nx) 
        idx -= y*nx
        x = idx
        vol[y,x,z] = vol_flat[id]
    return vol

if __name__ == "__main__":
    pass
