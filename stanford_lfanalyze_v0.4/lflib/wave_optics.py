import numpy as np

# --------------------------------------------------------------------------
#                           EMITTERS
# --------------------------------------------------------------------------

def spherical_wavefront(amplitude, wavelength, medium_index, L, M, x, y, z):
    '''
    Creates a 2D wavefront emanating from a point source at (x,y,z).
    This function does not make use of a paraxial approximation.
    '''
    dx = L/M
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    y_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    X, Y = np.meshgrid(x_ticks - x, y_ticks - y);

    r = np.sqrt(np.square(z) + np.square(X) + np.square(Y))
    k = medium_index * 2 * np.pi / wavelength
    u = amplitude / r * np.exp(-np.sign(z) * 1j * k * r)

    return u

# --------------------------------------------------------------------------
#                    AMPLITUDE & PHASE MASKS
# --------------------------------------------------------------------------

def circ_aperture(u_in, L, x, y, diameter):
    '''
    Apply an aperture mask of a given diameter to the wavefront.
    '''
    (M, N) = u_in.shape

    dx = L/M;
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    y_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    X, Y = np.meshgrid(x_ticks - x, y_ticks - y);
    R = np.sqrt(np.square(X) + np.square(Y));

    u_out = u_in;
    u_out[ np.nonzero( R >= (diameter / 2.0)) ] = 0;
    return u_out

def quadratic_phase(u_in, L, wavelength, medium_index, d_tubelens, f_tubelens):
    '''
    Apply an aperture mask of a given diameter to the wavefront.
    '''
    (M, N) = u_in.shape

    k = 2*np.pi*medium_index/wavelength;

    dx = L/M;
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    X, Y = np.meshgrid(x_ticks, x_ticks);

    magphase = np.exp( 1j * k / (2*f_tubelens) * (1-d_tubelens/f_tubelens)*( np.square(X)+np.square(Y) ));
    return u_in * magphase;

def abbe_sine_apodization(u_in, L, x, y, back_aperture_diameter, wavelength, medium_index, f_objective):
    '''
    Apply a the apodization function for the Abbe sine condition.

    Note: only use this at the Fourie plane!
    '''
    (M, N) = u_in.shape

    dx = L/M;
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    y_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    X, Y = np.meshgrid(x_ticks-x, y_ticks-y);
    R = np.sqrt(np.square(X) + np.square(Y));
    R[np.nonzero(R>(back_aperture_diameter/2.0))] = 0;

    thetas = np.arcsin(R*wavelength/medium_index);
    apodization = np.sqrt(np.cos(thetas));
    return apodization * u_in;


# --------------------------------------------------------------------------
#                              LENSES
# --------------------------------------------------------------------------

def ulens_array(uin,L,wavelength,D,ns,nu,zf,aperture_type='rect',fill_factor=1.0):
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
    aperture_type - the shape of the aperture (either 'rect' or 'circ')
    fill_factor - the fill factor of a ulens (from 0 to 1.0)
    uout - output field
    '''
    if aperture_type != 'rect' and aperture_type != 'circ':
        print 'ERROR: \'',aperture_type,'\' is an unrecognized microlens profile.  Choose \'rect\' or \'circ\''

    (M,N) = uin.shape
    dx=L/M;
    lenslet_len = int(M/ns);

    # Create a local optical system that retards the phase for each lenslet
    k = 2 * np.pi / wavelength;
    x = np.linspace(-D/2+dx/2, D/2-dx/2, nu)
    X, Y = np.meshgrid(x,x)
    phase = np.exp(-1j*k/(2*zf)*( np.square(X)+np.square(Y) )); # apply focus

    # Build the aperture mask - with ssf times the resolution of uin, blur and
    # sample back to uin resolution
    ssf=10
    nu_ss = nu*ssf
    dx_ss = dx/ssf
    x_ss = np.linspace(-D/2+dx_ss/2, D/2-dx_ss/2, nu_ss)
    b = np.ones((ssf,ssf)) / np.square(ssf)
    if (aperture_type == 'rect'):
        aperture_1d = np.logical_and(x_ss >= -D*fill_factor/2+dx_ss/2, x_ss <= D*fill_factor/2-dx_ss/2)
        aperture_2d = np.outer(aperture_1d, aperture_1d)
    elif (aperture_type == 'circ'):
        rad = D*fill_factor/2
        X_ss, Y_ss = np.meshgrid(x_ss,x_ss)
        R_ss = np.sqrt(np.square(X_ss) + np.square(Y_ss))
        aperture_2d = R_ss <= rad

    import scipy.ndimage as spndi
    aperture = spndi.convolve(aperture_2d.astype(np.float64),b,mode='nearest')
    hssf = np.floor(ssf/2)
    aperture = aperture[(hssf-1)::ssf, (hssf-1)::ssf];

    uout = uin.copy()
    for s in np.arange(ns):
        for t in np.arange(ns):
            cell = uin[s*lenslet_len:(s+1)*lenslet_len, t*lenslet_len:(t+1)*lenslet_len];
            if fill_factor < 1.0 or aperture_type == 'circ':
                uout[s*lenslet_len:(s+1)*lenslet_len, t*lenslet_len:(t+1)*lenslet_len] = cell * phase * aperture;
            else:
                uout[s*lenslet_len:(s+1)*lenslet_len, t*lenslet_len:(t+1)*lenslet_len] = cell * phase;
    return uout

# --------------------------------------------------------------------------
#                                PROPAGATORS
# --------------------------------------------------------------------------

def propTF(u1,L,wavelength,z):
    '''
    Propagate a wavefield a distance 'z' using the Fresnel transfer function method.

    This code was adapted from 'Computational Fourier Optics: a MATLAB Tutorial' by David Voelz
    
    propagation - transfer function approach
    assumes same x and y side lengths and
    uniform sampling
    u1 - source plane field
    L - source and observation plane side length
    wavelength - wavelength
    z - propagation distance
    u2 - observation plane field
    '''
    
    (M,N) = u1.shape

    dx=L/M;
    k = 2*np.pi/wavelength;

    # Compute the frequency sampling rate and frequency bin size.
    fs = 1.0/dx
    df = fs / M

    # The frequencies for an FFT are computed from the number of samples
    # and sampling rate, and it depends on whether the sequence length
    # is even or odd.
    #
    # See:
    #    
    #  http://www.astro.virginia.edu/class/oconnell/astr511/idl_5.1_html/idl9d.htm
    #  http://northstar-www.dartmouth.edu/doc/idl/html_6.2/FFT.html
    #
    if np.mod(M,2) == 0:       # Frequencies for FFT with even number of samples
        fx = np.concatenate((
            np.arange(0, (M/2.0-1.0) + 1) * df,
            np.array([fs/2]),
            np.arange(-(M/2.0-1.0), 0) * df
            ))
    else:                      # Frequencies for FFT with odd number of samples
        fx = np.concatenate((
            np.arange(0, (M/2.0-0.5) + 1) * df,
            np.arange(-(M/2.0-0.5), 0) * df
            ))

    FX,FY = np.meshgrid(fx,fx);
    H = np.exp(-1j*np.pi*wavelength*z*(np.square(FX)+np.square(FY)));  # Transfer function
  
    U1 = np.fft.fft2(u1);
    U2 = H * U1;
    u2 = np.fft.ifft2(U2);        # inv fft, center obs field

    return u2

def fftpad(u1):
    temp = np.zeros((u1.shape[0]*2, u1.shape[1]*2), dtype = np.complex)
    temp[0:u1.shape[0],0:u1.shape[1]] = u1; # Zero pad the input
    return temp

def fftunpad(u1):
    (M,N) = u1.shape
    return u1[0:M/2, 0:N/2];                 # remove zero padding


# ------------------------------------------------------------------------------------------------
#                           AMPLITUDE POINT SPREAD FUNCTION (APSF)
# ------------------------------------------------------------------------------------------------


def complex_quadrature(func, a, b, **kwargs):
    '''
    Allows us to use scipy numerical integration routines with a
    function that returns complex values.
    '''
    import scipy
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x):
        return scipy.imag(func(x))

    from scipy.integrate import quad
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])




def apsf_analytic(intensity, L, M, wavelength, objective_na, objective_mag, medium_index, f_tubelens, x, y, z):
    '''
    This code implements Equation 3.4.15 in 'Advanced Optical Imaging
    Theory' by Min Gu.  It gives the complex amplitude at the image plane of the
    microscope generated by an ideal point source at position x, y, z.

    Input parameters:

      intensity - intensity of the point source emitter
      L - linear dimension (used for width & height) in meters of the PSF simulation
      M - the number of pixel in each linear dimension
      lambda - wavelength of the simulated light emitted by the point source (meters)
      NA - numerical aperture of the objective
      mag - magnification of the objective
      medium_index - index of refraction of the sample medium
      f_tubelens - tube lens focal length (e.g. 200mm for Nikon, 180mm for Olympus, etc.)
      x,y,z - position of the point source

    Returns a tuple containing:

      apsf_img - 2D image of the APSF at the sensor behind the microlenses
      apsf_idxs - a list of radial coordinates
      apsf_valuse - radial values of the psf at the above coordinates
    
    '''
    
    dx = L/M;

    f_objective = f_tubelens/objective_mag;           # Objective focal length
    d1 = f_objective;                                 # d1 = f1 ?? (I'm not sure of this...)
    alpha_o = np.arcsin(objective_na/medium_index);   # Object side numerical aperture angle
    k = medium_index * 2 * np.pi / wavelength;        # Wave number
  
    # The Fresnel integral function is suitable for numerical apertures < ~0.5
    def fresnel_integral_fn(rho, u, v):
        import scipy.special
        return np.exp(1j*u/2*np.square(rho)) * scipy.special.jn(0, rho * v) * 2 * np.pi * rho;

    # This Debye integral (based on eqn. 6.2.17) is suitable for numerical apertures > ~0.5.
    #
    # We use a abbe sine apodization function here: P(theta) = P(r)*sqrt(cos(theta))
    def debye_integral_fn(theta, alpha, u, v):
        import scipy.special
        return np.sqrt(np.cos(theta))*np.exp(1j*u*np.square(np.sin(theta/2)) /
                                             (2 * np.square(np.sin(alpha/2))) ) * scipy.special.jn(0, np.sin(theta)/np.sin(alpha)*v)*np.sin(theta);

    def apsf_fn(u,v):
        if objective_na <= 0.4:                  # Fresnel (i.e. paraxial) theory for low NA objectives
            cmplx_integral = complex_quadrature(lambda rho: fresnel_integral_fn(rho, u, v),
                                                0.0, 1.0, limit=200)[0]
        else:                                     # Debye theory for high NA objectives
            cmplx_integral = complex_quadrature(lambda theta: debye_integral_fn(theta, alpha_o, u, v),
                                                0.0, 1.0, limit=200)[0]
        return np.exp( -1j*u / ( 4*np.square(np.sin(alpha_o/2)) ) ) * cmplx_integral;

        # These additional intensity scaling terms would seem not to
        # be necessary for us, though they appear in Gu's book.
        # Leaving them out for now... -broxton
        #
        #return (objective_mag / (np.square(d1) * np.square(wavelength)) *
        #        np.exp( -1j*u / ( 4*np.square(np.sin(alpha_o/2)) ) ) * cmplx_integral);

    # Compute range of radius values needed to compute the image
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, 2*M)
    X, Y = np.meshgrid(x_ticks - x, x_ticks - y);
    R = np.sqrt(np.square(X)+np.square(Y));
    min_r = R.min()
    max_r = R.max()
  
    # Create a lookup table of radial PSF samples.  Sample at twice the
    # grid frequency to make sure we can intepolate reasonably below.
    apsf_idxs = np.linspace(0, max_r+dx/2, M)

    # Compute the APSF for each element of the lookup table.
    apsf_vals = np.zeros(apsf_idxs.shape, dtype=np.complex);
    for i in range(len(apsf_idxs)):
        v = k * apsf_idxs[i] * np.sin(alpha_o);
        u = 4 * k *  z * np.square(np.sin(alpha_o/2));
        apsf_vals[i] = np.sqrt(intensity) * apsf_fn(u,v);  

    # Interpolate to populate the image.
    apsf_img = apsf_eval(L, M, apsf_idxs, apsf_vals, x, y);
    return (apsf_img, apsf_idxs, apsf_vals)


def apsf_fft(intensity, L, M, wavelength, objective_na,
             objective_mag, medium_index, f_tubelens, x, y, z):
    '''
    This code implements compute a microscope APSF using a Fourier
    transform method.  This is faster, though possibly a little bit
    less accurate than the method implemented in apsf_analytic.m

    Input parameters:

      intensity - intensity of the point source emitter
      L - linear dimension (used for width & height) in meters of the PSF simulation
      M - the number of pixel in each linear dimension
      lambda - wavelength of the simulated light emitted by the point source (meters)
      NA - numerical aperture of the objective
      mag - magnification of the objective
      medium_index - index of refraction of the sample medium
      f_tubelens - tube lens focal length (e.g. 200mm for Nikon, 180mm for Olympus, etc.)
      x,y,z - position of the point source

    Returns a tuple containing:

      apsf_img - 2D image of the APSF at the sensor behind the microlenses
      apsf_idxs - a list of radial coordinates
      apsf_valuse - radial values of the psf at the above coordinates
    
    '''
  
    # Handle points *on* the native plane by pushing them slightly off
    # of it by a sub-diffraction amount.  (The algorithm below does
    # not work if z = 0)
    if np.abs(z) < 0.1e-6:
      z = 0.1e-6;
  
    # Generate a spherical wavefront on the optical axis
    u1 = spherical_wavefront(np.sqrt(intensity), wavelength, medium_index, L, M, 0.0, 0.0, z); 

    # Apply the back aperture (i.e. NA limiting) function in the
    # fourier plane.  The FFT should be padded here since the point
    # source reaches out to the very edge of the simulation, and we
    # don't want any ringing as it wraps around circularly.
    U1 = np.fft.fftshift(np.fft.fft2(fftpad(u1)));                            # xform to aperture plane

    # Compute the sampling rate at the native plane, and then the size
    # and sampling rate and size of k-space.
    Lfourier = M / L;

    # Compute the back aperture diameter (in units of spatial frequencies - 1/m)
    back_aperture_diameter = 2 * objective_na / wavelength                     # From equation 21-4 in Gross

    # The center of the back aperture may not be exactly at the center of the FFT for even field sizes. 
    offset_x = 0;  offset_y = 0;
    if np.mod(U1.shape[0],2) == 0:
        offset_x = Lfourier/(4*M);
        offset_y = Lfourier/(4*M);

    # Apply the amplitude mask of the back aperture, and then apply an
    # apodization appropriate to account for an objective with the abbe
    # sine correction.
    U2 = circ_aperture(U1, Lfourier, offset_x, offset_y, back_aperture_diameter);   
    U2 = abbe_sine_apodization(U2, Lfourier, offset_x, offset_y, back_aperture_diameter, wavelength,
                               medium_index, f_tubelens/objective_mag);
    
    # Transform back to get the field at the image plane.
    u2 = fftunpad(np.fft.ifft2(np.fft.fftshift(U2)));                        # xform to image plane

    # Trim the edges off the FFT, since there are edge artifacts there.
    u2 = circ_aperture(u2, L, 0.0, 0.0, 0.98*L);   
  
    # Interpolate to populate the image.
    dx = L/M;
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    apsf_idxs = {}
    apsf_idxs['X'] = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    apsf_idxs['Y'] = np.linspace(-L/2+dx/2, L/2-dx/2, M)
    apsf_vals = u2;
    apsf_img = apsf_eval(L, M, apsf_idxs, apsf_vals, x, y);

    return (apsf_img, apsf_idxs, apsf_vals)



def apsf_eval(L, M, apsf_idxs, apsf_vals, x, y):
    '''
    This function takes an APSF lookup up table (computed using
    apsf_analytic.m or apsf_fft.m) and uses the translational invariance
    property of the APSF to interpolate and place it at any
    location x, y at the depth z for which apsf_idxs and apsf_vals was computed.
    '''
    
    dx = L/M;
    x_ticks = np.linspace(-L/2+dx/2, L/2-dx/2, M)

    # If this is a 2D lookup table, interpolate in 2D.
    if type(apsf_idxs) is dict:

        # Use 2D interpolation for FFT generate PSF data.
        apsf_idxs_x = apsf_idxs['X'];
        apsf_idxs_y = apsf_idxs['Y'];
        from scipy.ndimage import shift
        apsf_img_real = shift(np.real(apsf_vals), (float(y)/dx, float(x)/dx))
        apsf_img_imag = shift(np.imag(apsf_vals), (float(y)/dx, float(x)/dx))
        return apsf_img_real + 1j * apsf_img_imag;

    else:

        # Compute the range of
        X, Y = np.meshgrid(x_ticks - x, x_ticks - y)
        R = np.sqrt(np.square(X)+np.square(Y));

        # Interpolate to populate the image.
        from scipy.interpolate import griddata
        #    apsf_img =  griddata(apsf_idxs, apsf_vals, R, 'linear', 0.0);  # Fill value = 0.0
        apsf_img_real = griddata(apsf_idxs, np.real(apsf_vals), R, 'cubic', 0.0);  # Fill value = 0.0
        apsf_img_imag = griddata(apsf_idxs, np.imag(apsf_vals), R, 'cubic', 0.0);  # Fill value = 0.0
        return apsf_img_real + 1j * apsf_img_imag;


# ------------------------------------------------------------------------------------------------
#                                FULL LIGHT FIELD APSF
# ------------------------------------------------------------------------------------------------

def light_field_apsf(sim_size_m, sim_size_px, wavelength, objective_magnification,
                     objective_na, f_tubelens, medium_index, ulens_focal_length,
                     ulens_focal_distance, ulens_pitch,
                     ulens_profile, ulens_fill_factor,
                     ns, nu, x, y, apsf_idxs, apsf_vals):
    '''
    This routine computes a light field "splat" function using either
    Equation 3.4.15 in Advanced Optical Imaging Theory by Min Gu or
    using a faster FFT-based method.
    
    This method is slower than the FFT method below, but works well at
    all z-planes.
    '''

    # Compute defocus APSF of a point source at the image plane.
    u2 = apsf_eval(sim_size_m, sim_size_px, apsf_idxs, apsf_vals, x, y);

    # Phase curvature
    # u2 = quadratic_phase(u2, sim_size_m * objective_magnification, wavelength, 1.0, d_tubelens, f_tubelens);

    # Apply microlens phase mask
    u3 = ulens_array(u2, sim_size_m * objective_magnification, wavelength,
                     ulens_pitch, ns, sim_size_px / ns, ulens_focal_length,
                     aperture_type = ulens_profile, fill_factor = ulens_fill_factor);

    # Propagate from lenslets to sensor
    return propTF(u3, sim_size_m * objective_magnification, wavelength, ulens_focal_distance);


