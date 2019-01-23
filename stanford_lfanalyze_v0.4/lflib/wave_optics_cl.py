import sys
import numpy as np

try:
    from pyfft.cl import Plan
    import pyopencl as cl
    import pyopencl.array as cl_array
    LFLIB_HAVE_OPENCL = True
except ImportError:
    LFLIB_HAVE_OPENCL = False

def nextPowerOfTwo(value):
    ''' Determine how far past the nearest multiple of the value.
    This is useful for padding image dimensions to nice powers of 16,
    which helps promote efficient memory access on the GPU.'''

    exponent = np.log2(value)
    return int(np.power(2, np.ceil(exponent)));

src = '''
#include <pyopencl-complex.h>

__kernel void spherical_wavefront(__global float *X_mesh, __global float *Y_mesh,
                                  float A, float k, float x, float y, float z,
                                  __global cfloat_t *result) {

  // Grab the idx of the work item
  int i = get_global_id(0);

  // Create the spherical wavefront
  float xx = X_mesh[i]-x;
  float yy = Y_mesh[i]-y;
  float r = sqrt(z*z + xx*xx + yy*yy);
  float angle = -z / fabs(z) * k * r;
  result[i].x = A / r * cos(angle);  // Real part
  result[i].y = A / r * sin(angle);  // Imaginary part
}

__kernel void circ_aperture(__global float *X_mesh, __global float *Y_mesh, __global cfloat_t *U,
                            float diameter, float x, float y) {

  // Grab the idx of the work item
  int i = get_global_id(0);

  // Apply the aperture function
  float xx = X_mesh[i]-x;
  float yy = Y_mesh[i]-y;
  float r = sqrt(xx*xx + yy*yy);
  if (r >= diameter / 2.0) {
    U[i].x = 0;
    U[i].y = 0;
  } 
}

__kernel void ulens_array(__global float *X_mesh, __global float *Y_mesh, __global cfloat_t *U,
                          float k, float zf, int ns, int nu,
                          int width, int height) {
  // Grab the idx of the work item
  int x = get_global_id(0);        int y = get_global_id(1);
  int u = x % nu;                  int v = y % nu;

  float xx = X_mesh[v*nu+u];
  float yy = Y_mesh[v*nu+u];
  float angle = -k / (2 * zf) * (xx*xx + yy*yy);
  cfloat_t phase;
  phase.x = cos(angle);
  phase.y = sin(angle);

  // Grab the idx of the work item
  U[y*width+x] = cfloat_mul(U[y*width+x], phase);
}


__kernel void complex_multiply(__global cfloat_t *U1, __global cfloat_t *U2) {

  // Grab the idx of the work item
  int i = get_global_id(0);

  U1[i] = cfloat_mul(U1[i], U2[i]);
}

__kernel void intensity(__global cfloat_t *Uin, __global float *Iout) {

  // Grab the idx of the work item
  int i = get_global_id(0);

  cfloat_t a = Uin[i];
  Iout[i] = (a.x*a.x) - (a.y*a.y);
}

__kernel void fftpad(__global cfloat_t *unpadded_buf,
                     __global cfloat_t *padded_buf,
                     int unpadded_x_stride, int padded_x_stride) {

  // Grab the idx of the work item
  int x = get_global_id(0);
  int y = get_global_id(1);
  int unpadded_idx = y * unpadded_x_stride + x;
  int padded_idx = y * padded_x_stride + x;

  padded_buf[padded_idx] = unpadded_buf[unpadded_idx];
}

__kernel void fftunpad(__global cfloat_t *padded_buf,
                       __global cfloat_t *unpadded_buf,
                       int padded_x_stride, int unpadded_x_stride) {

  // Grab the idx of the work item
  int x = get_global_id(0);
  int y = get_global_id(1);
  int padded_idx = y * padded_x_stride + x;
  int unpadded_idx = y * unpadded_x_stride + x;

  unpadded_buf[unpadded_idx] = padded_buf[padded_idx];
}

__kernel void fftshift(__global cfloat_t *buf,
                       __global cfloat_t *shifted_buf,
                       int width, int height) {

  int p = ceil(width / 2.0);
  int q = ceil(height / 2.0);

  // Grab the idx of the work item
  int x = get_global_id(0);
  int y = get_global_id(1);
  int idx = y * width + x;
  int shift_idx;
  if ((y >= q) && (x >= p))
    shift_idx = (y-q) * width + (x-p);

  else if ((y >= q) && (x < p))
    shift_idx = (y-q) * width + (x+p);

  else if ((y < q) && (x >= p))
    shift_idx = (y+q) * width + (x-p);

  else 
    shift_idx = (y+q) * width + (x+p);

  shifted_buf[shift_idx] = buf[idx];
}


'''

cl_platform = cl.get_platforms()[0]

class WaveOpticsCl(object):
    def __init__(self, gpu_id = 5):

        if LFLIB_HAVE_OPENCL:
            # Set up OpenCL
            from pyopencl.tools import get_gl_sharing_context_properties

            print 'init started'            
            if sys.platform == "darwin":
                self.cl_ctx = cl.Context() 
            elif gpu_id is not None:
                self.cl_ctx = cl.Context(properties=[(cl.context_properties.PLATFORM, cl_platform)],
                                         devices = [cl_platform.get_devices()[gpu_id]])
            else:
                self.cl_ctx = cl.Context(properties=[(cl.context_properties.PLATFORM, cl_platform)])
            print 'compiling and queue generating'
            self.cl_prog = cl.Program(self.cl_ctx, src).build()
            self.cl_queue = cl.CommandQueue(self.cl_ctx)
            self.pyfft_plan = None
            self.plan_shape = None
            self.plan_dtype = None
            print 'init_ended'
            
        else:  # No OpenCL
            print "ERROR: your platform does not seem to support OpenCL.  These are required to generate wavespreads.   Exiting."
            raise SystemExit


    # --------------------- FOURIER TRANSFORM ROUTIENS (PYFFT) -----------------------

    def fft(self, gpu_data):

        # Create a new pyfft plan, if necessary
        if self.plan_shape != gpu_data.shape or self.plan_dtype != gpu_data.dtype:
            self.pyfft_plan = Plan(gpu_data.shape, dtype = gpu_data.dtype, queue=self.cl_queue, fast_math = False)
            self.plan_dtype = gpu_data.dtype
            self.plan_shape = gpu_data.shape

        self.pyfft_plan.execute(gpu_data.data)
        return gpu_data

    def ifft(self, gpu_data):

        # Create a new pyfft plan, if necessary
        if self.plan_shape != gpu_data.shape or self.plan_dtype != gpu_data.dtype:
            self.pyfft_plan = Plan(gpu_data.shape, dtype = gpu_data.dtype, queue=self.cl_queue, fast_math = False)
            self.plan_dtype = gpu_data.dtype
            self.plan_shape = gpu_data.shape

        self.pyfft_plan.execute(gpu_data.data, inverse = True)
        return gpu_data

    # --------------------------- WAVE OPTICS ROUTINES ---------------------------------

    def prepare_meshgrid(self, L, M):
        dx = L/M
        x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
        y_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
        X, Y = np.meshgrid(x_ticks, y_ticks);
        X_gpu = cl_array.to_device(self.cl_queue, X.astype(np.float32))
        Y_gpu = cl_array.to_device(self.cl_queue, Y.astype(np.float32))
        return (X_gpu, Y_gpu)

    def fftpad(self, unpadded_buf, pad2x = True):
        if pad2x:
            padded_shape = (nextPowerOfTwo(unpadded_buf.shape[0]*2), nextPowerOfTwo(unpadded_buf.shape[1]*2))
        else:
            padded_shape = (nextPowerOfTwo(unpadded_buf.shape[0]), nextPowerOfTwo(unpadded_buf.shape[1]))
        padded_buf = cl_array.zeros(self.cl_queue, padded_shape, dtype = unpadded_buf.dtype)
        self.cl_prog.fftpad(self.cl_queue, unpadded_buf.shape, None,
                            unpadded_buf.data, padded_buf.data,
                            np.uint32(unpadded_buf.shape[1]),
                            np.uint32(padded_buf.shape[1]))
        return padded_buf

    def fftunpad(self, padded_buf, unpadded_shape):
        unpadded_buf = cl_array.zeros(self.cl_queue,
                                      unpadded_shape,
                                      dtype = padded_buf.dtype)
        self.cl_prog.fftunpad(self.cl_queue, unpadded_buf.shape, None,
                              padded_buf.data, unpadded_buf.data,
                              np.uint32(padded_buf.shape[1]),
                              np.uint32(unpadded_buf.shape[1]))
        return unpadded_buf

    def fftshift(self, buf):
        shifted_buf = cl_array.empty_like(buf)
        self.cl_prog.fftshift(self.cl_queue, buf.shape, None,
                              buf.data, shifted_buf.data,
                              np.uint32(buf.shape[1]), np.uint32(buf.shape[0]))
        return shifted_buf

    def spherical_wavefront(self, X_mesh, Y_mesh, intensity, wavelength, x, y, z, result_buf = None):
        '''
        Creates a 2D wavefront emanating from a point source at (x,y,z).
        This function does not make use of a paraxial approximation...

        Supply an existing gpu buffer as result_buf to make things
        slightly faster.
        '''
        if result_buf == None:
            result_buf = cl_array.empty(self.cl_queue, X_mesh.shape, dtype=np.complex64)

        k = 2 * np.pi / wavelength
        num_pixels = int(np.prod(X_mesh.shape))
        self.cl_prog.spherical_wavefront(self.cl_queue, (num_pixels,), None,
                                         X_mesh.data, Y_mesh.data,
                                         np.float32(intensity), np.float32(k),
                                         np.float32(x), np.float32(y), np.float32(z), result_buf.data)
        return result_buf

    def circ_aperture(self, X_mesh, Y_mesh, u_in, x, y, diameter):
        '''
        Apply an aperture mask of a given diameter to the wavefront.

        NOTE: changes the wavefront in place.
        '''
        num_pixels = int(np.prod(X_mesh.shape))
        self.cl_prog.circ_aperture(self.cl_queue, (num_pixels,), None,
                                   X_mesh.data, Y_mesh.data, u_in.data,
                                   np.float32(diameter), np.float32(x), np.float32(y))
        return u_in

    def ulens_array(self, X_mesh, Y_mesh, u_in, wavelength, ns, nu, zf):
        '''
        Apply a tiled phase mask created by a lenslet array.  Be sure that
        the wavefront you pass in has dimensions that are equal to ns *
        nu.

        uin - input field
        L - side length
        wavelength - wavelength
        D - lenslet pitch
        ns - number of lenslets in each linear dimension
        zf - focal distance (+ converge, - diverge)
        uout - output field

        NOTE: changes the wavefront in place.
        '''
        k = 2 * np.pi / wavelength
        self.cl_prog.ulens_array(self.cl_queue, u_in.shape, None,
                                 X_mesh.data, Y_mesh.data, u_in.data,
                                 np.float32(k), np.float32(zf), np.int32(ns), np.int32(nu),
                                 np.int32(u_in.shape[1]), np.int32(u_in.shape[0]))
        return u_in


    def propTF(self, u1_buf, L ,wavelength, z):
        '''
        propagation - transfer function approach
        assumes same x and y side lengths and
        uniform sampling
        u1 - source plane field
        L - source and observation plane side length
        wavelength - wavelength
        z - propagation distance
        u2 - observation plane field
        '''
        # Create a padded U1 buffer
        u1_padded_buf = self.fftpad(u1_buf, pad2x = False)
        Lpad = L * u1_padded_buf.shape[0] / u1_buf.shape[0]
    
        (M,N) = u1_padded_buf.shape
        dx=Lpad/M;
        k=2*np.pi/wavelength;

        # Compute the proper frequency coordinates for even vs. odd field size
        fs = 1.0/dx
        df = fs / M
        fx = np.linspace(0,fs-df,M) - (fs - np.mod(M,2)*df) / 2
        FX,FY = np.meshgrid(fx,fx);
        H = np.exp(-1j*np.pi*wavelength*z*(np.square(FX)+np.square(FY)));  # Transfer function
        H = np.fft.fftshift(H);

        # Send the transfer function to the GPU
        H_buf = cl_array.to_device(self.cl_queue, H.astype(np.complex64))

        U1_buf = self.fft(self.fftshift(u1_padded_buf))
        num_pixels = int(np.prod(U1_buf.shape))
        self.cl_prog.complex_multiply(self.cl_queue, (num_pixels,), None,
                                      U1_buf.data, H_buf.data)
        u2_buf = self.ifft(U1_buf)

        result = u2_buf.get(self.cl_queue).reshape(u2_buf.shape)
        u2 = np.fft.ifftshift(result);                   # inv fft shift
        return u2[0:u1_buf.shape[0], 0:u1_buf.shape[1]]  # Remove padding

