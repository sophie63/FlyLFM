# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import os, sys

import array
import numpy as np

from lflib.imageio import load_image, save_image
from lflib.lightfield import LightField


#----------------------------------------------------------------------------------
#                         RAYSPREAD SUPPORT ROUTINES
#----------------------------------------------------------------------------------

class RayspreadDatabase(object):

    def __init__(self, nu, nv, ns, nt, z_center, z_coords, supersample_factor):
        self.rayspreads = {}
        self.flipped_rayspreads = {}
        self.nu = nu
        self.nv = nv
        self.ns = ns
        self.nt = nt
        self.nrays = self.nu*self.nv*self.ns*self.nt
        self.z_coords = z_coords
        self.nx = ns * supersample_factor
        self.ny = nt * supersample_factor
        self.nz = len(z_coords)
        self.nvoxels = self.nx*self.ny*self.nz
        self.supersample_factor = supersample_factor
        self.z_center = z_center

    def empty_volume(self):
        '''Returns an empty volume with the correct dimensions to match this raydb.'''
        return  np.zeros((self.ny, self.nx, self.nz), dtype=np.float32)

    def empty_lf(self):
        '''Returns an empty lightfield with the correct dimensions to match this raydb.'''
        lf = np.zeros((self.nv*self.nt, self.nu*self.ns), dtype=np.float32)
        return LightField(lf, self.nu, self.nv, self.ns, self.nt, representation = LightField.TILED_LENSLET)

    def ones_volume(self):
        '''Returns an empty volume with the correct dimensions to match this raydb.'''
        return  np.ones((self.ny, self.nx, self.nz), dtype=np.float32)

    def ones_lf(self):
        '''Returns an empty lightfield with the correct dimensions to match this raydb.'''
        lf = np.ones((self.nv*self.nt, self.nu*self.ns), dtype=np.float32)
        return LightField(lf, self.nu, self.nv, self.ns, self.nt, representation = LightField.TILED_LENSLET)

def compute_rayspread_imap(in_tuple):
    '''
    This is a helper routine that allows us to build the rayspread
    database in a multithreaded fashion.
    '''
    from pylflib import compute_rayspread
    spread = compute_rayspread(in_tuple[0], in_tuple[1], in_tuple[2], in_tuple[3], in_tuple[4], in_tuple[5],
                               sine_correction = True,
                               supersample_factor = in_tuple[6], num_raysamples = in_tuple[8],
                               num_zsamples = in_tuple[9],
                               pixel_fill_factor = in_tuple[10], ulens_fill_factor = in_tuple[11],
                               ulens_profile = in_tuple[12])
    return (in_tuple, spread)

#----------------------------------------------------------------------------------
#                         WAVESPREAD SUPPORT ROUTINES
#----------------------------------------------------------------------------------

class ApsfCache(object):
    '''
    A simple thread-safe dictionary object for storing APSFs.
    '''

    def __init__(self):
        from multiprocessing import Pool, Manager
        manager = Manager()                 # This is a strange work-around so that I can have a Lock()
        self.lock = manager.Lock()          # shared by a multi-processing pool...
        self.d = manager.dict()             # shared by a multi-processing pool...

    def has_key(self, key):
        self.lock.acquire()
        result = self.d.has_key(key)
        self.lock.release()
        return result

    def set(self, key, value):
        self.lock.acquire()
        self.d[key] = value
        self.lock.release()

    def get(self, key):
        self.lock.acquire()
        result = self.d[key]
        self.lock.release()
        return result

class PsfDatabase(object):

    def __init__(self, nu, nv, ns, nt, z_center, z_coords, supersample_factor):
        self.raw_splats = {}
        self.psf_coordinates = {}
        self.psf_coefficients = {}
        self.nu = nu
        self.nv = nv
        self.ns = ns
        self.nt = nt
        self.nrays = self.nu*self.nv*self.ns*self.nt
        self.z_coords = z_coords
        self.nx = ns * supersample_factor
        self.ny = nt * supersample_factor
        self.nz = len(z_coords)
        self.nvoxels = self.nx*self.ny*self.nz
        self.supersample_factor = supersample_factor
        self.z_center = z_center
        self.max_coefficient = 0.0;

    def empty_volume(self):
        '''Returns an empty volume with the correct dimensions to match this raydb.'''
        return  np.zeros((self.ny, self.nx, self.nz), dtype=np.float32)

    def empty_lf(self):
        '''Returns an empty lightfield with the correct dimensions to match this raydb.'''
        lf = np.zeros((self.nv*self.nt, self.nu*self.ns), dtype=np.float32)
        return LightField(lf, self.nu, self.nv, self.ns, self.nt, representation = LightField.TILED_LENSLET)

    def ones_volume(self):
        '''Returns an empty volume with the correct dimensions to match this raydb.'''
        return  np.ones((self.ny, self.nx, self.nz), dtype=np.float32)

    def ones_lf(self):
        '''Returns an empty lightfield with the correct dimensions to match this raydb.'''
        lf = np.ones((self.nv*self.nt, self.nu*self.ns), dtype=np.float32)
        return LightField(lf, self.nu, self.nv, self.ns, self.nt, representation = LightField.TILED_LENSLET)

def hann_window(x,y,z):
    xr = 0.0; yr = 0.0; zr = 0.0;

    xr = (np.cos(np.pi*x) + 1.0) * (np.abs(x) < 1.0);
    yr = (np.cos(np.pi*y) + 1.0) * (np.abs(y) < 1.0);
    zr = (np.cos(np.pi*z) + 1.0) * (np.abs(z) < 1.0);

    return xr * yr * zr;

def lanczos1_window(x,y,z):
    r = np.sqrt(x*x+y*y+z*z)
    if r > 1.0:
        weight = 0.0;
    else:
        weight = np.sinc(r) * np.sinc(r)
    return weight

def lanczos3_window(x,y,z):
    r = np.sqrt(x*x+y*y+z*z)
    if r > 3.0:
        weight = 0.0;
    else:
        weight = np.sinc(r) * np.sinc(r/3.0)
    return weight

def compute_light_field_psf(intensity, x, y, z,
                            z_spacing, num_lenslets_in_psf, supersample_factor,
                            lenslet_array, apsf_cache, scrub = False, compute_apsf_with_fft = True):

    '''
    Generates a light field "splat" function for a coherent point source at a given position (x,y,z).

    Parameters:

    intensity -- radiant intensity of the point source
    
    x,y,z     -- position of the point source

    num_lenslets_in_psf -- the number of lenslets to use when
                           simulating the PSF.  Be sure to specify a
                           large enough number of lenslets to contain
                           the entire splat function!

    lenslet_array       -- LensletArray object containing essential optical metadata

    wavelength_nm       -- Wavelength of light to be simulated

    scrub -- set to True to compute the incoherent sum of coherent PSFs over the area of
             of a voxel.  This is a useful post-aliasing filter when sampling the volume
             below the diffraction limit.

    compute_apsf_with_fft -- set to False to use the analytical FFT (much slower!)
    '''
    import time
    tic = time.time()
    
    # Extract optical parameters
    wavelength = lenslet_array.center_wavelength * 1e-9;
    objective_magnification = lenslet_array.objective_magnification
    objective_na = lenslet_array.objective_na
    medium_index = lenslet_array.medium_index
    ulens_pitch = lenslet_array.ulens_pitch * 1e-6
    ulens_focal_length = lenslet_array.ulens_focal_length * 1e-6
    ulens_focal_distance = lenslet_array.ulens_focal_distance * 1e-6
    f_tubelens = lenslet_array.tubelens_focal_length * 1e-3 # convert from mm to meters
        
    # Effective pitch of the microlens array when conjugated to the object side of the objective.
    effective_pitch = ulens_pitch / objective_magnification;  

    # SIMULATION PARAMS
    ns = num_lenslets_in_psf
    nu = lenslet_array.nu
    sim_size_m = effective_pitch * ns;   # Simulation size (meters)

    # We must be sure to run the simulation at or above the nyquist
    # rate.  The sampling rate is scaled here.
    sampling_rate = ulens_pitch / float(nu);
    nyquist_rate = wavelength * ulens_focal_length / (2.0 * ulens_pitch);
    sim_sampling = 2**int(np.floor(sampling_rate / nyquist_rate))
    # print '\t--> Upsampling the simulation by ', sim_sampling, 'x to reach nyquist rate.'
    
    sim_size_px = sim_sampling * ns * nu;               # Sim size in pixels

    # Compute "scrubbing" parameters
    #
    # Subsample each voxel.  For voxels near the native plane, we need
    # more subsamples in order to produce accurate wave optics psf
    # estimates.  Farther from the native plane, we can get away with
    # fewer subsamples.
    diffraction_limit = wavelength / (2 * lenslet_array.objective_na);

    # Scaling factor from normalized subsample coords to real coords in object space.
    k = effective_pitch / supersample_factor

    if scrub:
        if (ns <= 5):
            num_xy_samples = np.round(effective_pitch / (diffraction_limit * supersample_factor) ) * 4
            num_z_samples = 1
        elif (ns <= 9):
            num_xy_samples = np.round(effective_pitch / (diffraction_limit * supersample_factor) ) * 2
            num_z_samples = 1
        elif (ns <= 15):
            num_xy_samples = np.round(effective_pitch / (diffraction_limit * supersample_factor) )
            num_z_samples = 1
        else:
            num_xy_samples = np.round(effective_pitch / (diffraction_limit * supersample_factor) ) / 2
            num_z_samples = 1

        # Ensure we have at least one sample, and make sure the num_xy_samples is odd.
        if num_xy_samples < 3:
            num_xy_samples = 3

        if num_xy_samples % 2 == 0:
            num_xy_samples += 1

        x_subsamples = np.linspace( -1.0+1.0/num_xy_samples, 1.0-1.0/num_xy_samples, num_xy_samples )
        y_subsamples = np.linspace( -1.0+1.0/num_xy_samples, 1.0-1.0/num_xy_samples, num_xy_samples )
        if num_z_samples == 1:
            z_subsamples = np.array([0])
        else:
            z_subsamples = np.linspace( -1.0, 1.0, num_z_samples )

    else:  # scrubbing disabled
        x_subsamples = np.array([0])
        y_subsamples = np.array([0])
        z_subsamples = np.array([0])

    weight_sum = 0
    psf = None


    from wave_optics import apsf_analytic, apsf_fft, light_field_apsf
    for sz in z_subsamples:
        kz = z_spacing * 1e-6

        # Check to see if this APSF lookup is in our cache.  If not, compute it!
        if apsf_cache.has_key(z+kz*sz) == False:
            if compute_apsf_with_fft:
                # Intensity of point source (chosen to produce similar intensities to analytical method)
                intensity = 1e-12; 
                (apsf_img, apsf_idxs, apsf_vals) = apsf_fft(intensity, sim_size_m, sim_size_px, wavelength,
                                                            objective_na, objective_magnification,
                                                            medium_index, f_tubelens, 0, 0, z+kz*sz)
            else:
                # Intensity of point source (chosen to produce similar intensities to fft method)
                intensity = 35e2; 
                (apsf_img, apsf_idxs, apsf_vals) = apsf_analytic(intensity, sim_size_m, sim_size_px, wavelength,
                                                                 objective_na, objective_magnification,
                                                                 medium_index, f_tubelens, 0, 0, z+kz*sz)

            apsf_cache.set(z+kz*sz, (apsf_idxs, apsf_vals))

        # Extract the data from the cache, which should now be populated.
        (apsf_idxs, apsf_vals) = apsf_cache.get(z+kz*sz)

        # Give the voxel some aereal extent by "scrubbing" it around. (This avoids aliasing in the volume.)
        for sx in x_subsamples:
            for sy in y_subsamples:

                # Choose a weighting function manually.  Options here include: {lanczos1, hann, uniform}
                # weight = lanczos1_window( sx, sy, sz )                           # Lanczos-1
                weight = hann_window( sx, sy, sz )                                 # Hann
                # weight = (1.0-np.abs(sx)) * (1.0-np.abs(sy)) * (1.0-np.abs(sz))  # Bilinear

                U_sensor = weight * light_field_apsf(sim_size_m, sim_size_px, wavelength,
                                                     objective_magnification, objective_na, f_tubelens,
                                                     medium_index, ulens_focal_length, ulens_focal_distance,
                                                     ulens_pitch,
                                                     lenslet_array.ulens_profile, lenslet_array.ulens_fill_factor,
                                                     ns, nu, x+k*sx, y+k*sy, apsf_idxs, apsf_vals)
                new_psf = np.real(U_sensor*U_sensor.conj()).astype(np.float32)

                # Subsample to the true sensor resolution.  We
                # downsample by averaging blocks, since these blocks
                # are all recorded by the same square sensor pixel.
                # This actually produces some aliasing, but this
                # aliasing should be present in the real optical
                # system as well, so I think this is a physically
                # accurate averaging scheme. -broxton
                height, width = new_psf.shape
                new_psf = np.average(np.split(np.average(np.split(new_psf, width // sim_sampling, axis=1), axis=-1),
                                              height // sim_sampling, axis=1), axis=-1)
                
                
                weight_sum += weight

                # Add psf intensities incoherently!
                if psf != None:
                    psf += new_psf
                else:
                    psf = new_psf

    # Normalize by the weight sum and return the psf
    return psf / weight_sum

def compute_psf_imap(in_tuple):
    '''
    This is a helper routine that allows us to build the psf database in a multithreaded fashion.
    '''

    # Extract some basic optical parameters to use in the subsample interval calculations below.
    x = in_tuple[0]
    y = in_tuple[1]
    z = in_tuple[2]
    z_spacing = in_tuple[3]
    supersample_factor = in_tuple[4]
    num_lenslets_in_aperture = in_tuple[5]
    lenslet_array = in_tuple[6]
    radiometry = in_tuple[7]
    x_idx = in_tuple[8]
    y_idx = in_tuple[9]
    z_idx = in_tuple[10]
    voxels_as_points = in_tuple[11]
    apsf_cache = in_tuple[12]

    # PART 1:
    #
    # Compute light field PSFs.
    #

    # Turn on "scrubbing" which gives each voxel an areal extent via a simple numerical integral
    scrub = not voxels_as_points
    intensity = 1.0  # Arbitrary for now...
    psf = compute_light_field_psf(intensity, x, y, z, z_spacing, num_lenslets_in_aperture,
                                  supersample_factor, lenslet_array, apsf_cache, scrub)


    # The number of nonzero entries in a psf image is usually
    # small compared to the total number of pixels in that
    # image, so we extract those coordinates and store/process
    # them in a sparse representation.
    #
    nonzero_entries = np.nonzero(psf)
    coords = np.vstack((nonzero_entries[0], nonzero_entries[1])).T

    # Convert the coordinates from lenslet image coordinates
    # [m,n] to light field coordinates [u,v,s,t] and then
    # "recenter" the [s,t] coordinates so that the center
    # lenslet in the psf has coordinates [0,0].  This
    # convention allows us to know that the point source was
    # centered under lenslet [0,0], and this assumption will
    # be used in the forward_project() and back_project()
    # operators in volume.py
    #
    from lightfield import lenslet_image_coords_to_uvst
    uvst = lenslet_image_coords_to_uvst(coords, lenslet_array.nu, lenslet_array.nv)
    uvst[:,2] -= int((num_lenslets_in_aperture-1) / 2)
    uvst[:,3] -= int((num_lenslets_in_aperture-1) / 2)

    # Find the maximum coefficient across all PSFs.  We will
    # use this to trim out coefficients that would be
    # "undetectable" (i.e. below the quantization threshold of
    # our sensor) below.
    coefficients = psf[nonzero_entries]

    # PART 2:
    #
    # Eliminate "undetectable" coefficients that would fall below
    # the quantization level of our sensor.  Eliminating these
    # speeds up the forward/backward projection steps
    # considerably, and should not degrade reconstruction
    # performance.
    #

    # This performance shortcut avoids processing the rays
    # that are outside of the NA of the objective.  These rays
    # are vignetting and highly aberrated, and aren't so good
    # to use in any case.
    (uvst, coefficients) = lenslet_array.trim_to_na(uvst, coefficients)

    # Eliminate any coefficients that fall below the min detectable level for a 14 bit sensor.
    orig_intensity = coefficients.sum()
    num_nonzero = coefficients.shape[0]
    rejection_threshold = coefficients.sum() * 0.05        # Keep 95% of total intensity

    # Sort by increasing intensity
    keys = np.argsort(coefficients)
    uvst = uvst[keys,:]
    coefficients = coefficients[keys]

    # Compute the comulative sum, extract the index of
    # the first entry to keep.  Eliminate undetectable
    # entries from both coefficients and uvst matrices
    coef_cumsum = np.cumsum(coefficients)
    detectable_entries = np.nonzero(coef_cumsum > rejection_threshold)[0][0]
    coefficients = coefficients[detectable_entries:]
    uvst = uvst[detectable_entries:,:]

    # Sorting the coefficients so that identical [u,v]
    # pairs are grouped together will help increase
    # cache coherence when accessing the light field
    # texture.  This appears to speed things up in
    # deconvolution on the GPU by ~33%.
    #
    # Sort by s entries
    keys = np.argsort(uvst[:,2])
    uvst = uvst[keys,:]
    coefficients = coefficients[keys]

    # Sort by t entries (maintaining s order with merge sort)
    keys = np.argsort(uvst[:,3], kind='mergesort')
    uvst = uvst[keys,:]
    coefficients = coefficients[keys]

    # Sort by u entries
    keys = np.argsort(uvst[:,0], kind='mergesort')
    uvst = uvst[keys,:]
    coefficients = coefficients[keys]

    # Sort by v entries
    keys = np.argsort(uvst[:,1], kind='mergesort')
    uvst = uvst[keys,:]
    coefficients = coefficients[keys]

    pct_pixels_retained = 100.0 * coefficients.shape[0] / num_nonzero
    pct_intensity_retained = 100.0 * coefficients.sum() / orig_intensity

    # For debugging:
    #
    #psf_filename = 'psf_%d__%d_%d.tif' % (z_idx,x_idx,y_idx)
    #print 'SAVING ', psf_filename
    # invpsf_filename = 'invpsf_%d__%d_%d.tif' % (z,x,y)
    #from lflib.imageio import save_image
    #save_image(psf_filename, psf)
    # from lflib.lightfield import lenslet2subaperture
    # save_image(invpsf_filename,
    #            lenslet2subaperture(psf, lenslet_array.nu, lenslet_array.nv,
    #                                num_lenslets_in_aperture, num_lenslets_in_aperture))


    return (uvst, coefficients, x_idx, y_idx, z_idx, num_lenslets_in_aperture,
            psf.sum(), num_nonzero, pct_pixels_retained, pct_intensity_retained)


#----------------------------------------------------------------------------------
#                         LENSLET ARRAY CLASS                            
#----------------------------------------------------------------------------------

class LfSensor(object):
    def __init__(self, nu, nv, ns, nt):
        self.nu = nu; self.nv = nv;
        self.ns = ns; self.nt = nt;

    def uv_to_angle(self, u, v, sine_correction = True):
        """Convert from discrete ray angle direction to an angle."""
        raise NotImplementedError


class LensletArray(LfSensor):

    def __init__(self, nu, nv, ns, nt,
                 ulens_pitch, pixel_size, ulens_focal_length, ulens_focal_distance,
                 objective_magnification, objective_na, medium_index,
                 tubelens_focal_length,
                 ulens_fill_factor = 1.0, pixel_fill_factor = 1.0,
                 ulens_profile = 'rect',
                 center_wavelength = 509):  # Units: nanometers
        """
        Create a new lenslet array object.
        
        nu, nv, ns, nt - are the dimensions (in pixels) of the light
                          field (u,v for ray angles, s,t for spatial
                          extent)
        pitch          - pitch of the microlens array in microns.
                         (each microlens is ulens_pitch x ulens_pitch
                          micron in size)
        magnification  - magnification of the objective
        na             - numerical aperture of the objective
        medium_index   - index of refraction in the medium
        """
        LfSensor.__init__(self, nu, nv, ns, nt)
        self.ulens_pitch = ulens_pitch
        self.pixel_size = pixel_size
        self.ulens_focal_length = ulens_focal_length
        self.ulens_focal_distance = ulens_focal_distance
        self.objective_magnification = objective_magnification
        self.objective_na = objective_na
        self.medium_index = medium_index
        self.tubelens_focal_length = tubelens_focal_length
        self.ulens_fill_factor = ulens_fill_factor
        self.pixel_fill_factor = pixel_fill_factor
        self.ulens_profile = ulens_profile
        self.center_wavelength = center_wavelength
    
    def uv_to_angle(self, u, v, sine_correction = True):
        """
        Convert from discrete ray angle direction to an angle.  The
        lenslet class implements this use the Abbe sine correction.
        Returns a tuple (theta, phi) of angles relative to the optical
        axis of the lenslet.
        """
        
        if sine_correction:
            theta = -np.arcsin( 2 * (float(u)-(self.nu-1)/2.0) / (self.nu-1) * self.objective_na / self.medium_index )
            phi   = -np.arcsin( 2 * (float(v)-(self.nv-1)/2.0) / (self.nv-1) * self.objective_na / self.medium_index )

        else:
            theta = -2 * (float(u)-(self.nu-1)/2.0) / (self.nu-1) * np.arcsin( self.objective_na / self.medium_index )
            phi   = -2 * (float(v)-(self.nv-1)/2.0) / (self.nv-1) * np.arcsin( self.objective_na / self.medium_index )

        return (theta, phi)

    # --------------------
    # RAYSPREAD GENERATION
    # --------------------

    def rayspread_database(self, z_len, z_center, z_spacing, radiometry=None,
                           supersample_factor = 1, num_threads = 1, optimize = True):
        """
        Generate a rayspread database for this lenslet array based on
        the geometric optics model.  This database in needed by the
        reconstruction and projection algorithms in volume.py.

        z_len             - the number of slices in the stack
        z_center          - the z offset in micron from the focal plane
                            of the center slice in the volume
                            positive z is "up" from the focal plane towards
                            the objective.  If there are an even number of
                            slices in the volume, then this is the average
                            of the z offset of the two middle slices.
        z_spacing         - the number of micron between slices in the
                            volume
        radiometry        - is either 'None' or a mapping (u,v)->float 
                            that specifies the intensity of a uniform
                            fluorescent volume
                            along each of the ray angles.  All the (u,v)
                            ray angles that will be used in the reconstruction
                            must be specified in this mapping.
       supersample_factor - Request supersampled ray spread kernels.  This is
                            useful in super-resolved deconvolution.
        """
        # modify these to adjust the fill factor used by the rayspread code
        pixel_fill_factor = self.pixel_fill_factor
        ulens_fill_factor = self.ulens_fill_factor
        ulens_profile = self.ulens_profile

        # ---------------------
        # Compute z coordinates
        # ---------------------
        z_start = z_center - 0.5 * (z_len - 1) * z_spacing
        z_coords = [x * z_spacing + z_start for x in range(z_len)]

        # ---------------------
        # Compute Radiometry
        # ---------------------

        # If the user has not supplied radiometry, we compute "approximate radiometry here"
        if radiometry is None:
            print 'WARNING: Using approximate radiometry.  Supplying calibrated radiometry will yield better results.'
            radiometry = np.zeros((self.nv, self.nu), dtype=np.float32)
            dummy_radiometry = np.ones((self.nv, self.nu), dtype=np.float32)
            for v in range(self.nv):
                for u in range(self.nu):
                    #optimize = False
                    if not optimize:
                        # default sampling values                    
                        num_raysamples = 10
                        num_zsamples = 5
                    else:
                        # optimize raypspread computation by having a z-depth and angle dependent sampling rate
                        # (see Evernote 'rayspread code optimization')
                        #(theta, phi) = self.uv_to_angle(u, v)
                        z = 0
                        DX = self.ulens_pitch * 1.0 / self.objective_magnification / supersample_factor;
                        MIN_ZSAMPLE_FACTOR = 1.0; # increase this for finer z-sampling
                        THETA_NA = np.arcsin(self.objective_na / self.medium_index);
                        num_zsamples = np.ceil(z_spacing * MIN_ZSAMPLE_FACTOR / DX * np.tan(THETA_NA)); # this was set to "5" before
                        (theta1, phi1) =   self.uv_to_angle( u-pixel_fill_factor/2, v-pixel_fill_factor/2) # this is only approximate, since this uses the old incorrect uv_to_angle()
                        (theta2, phi2) =   self.uv_to_angle( u+pixel_fill_factor/2, v+pixel_fill_factor/2)
                        MIN_SAMPLE_FACTOR = 1.0; # increase this for finer (u,v) sampling
                        num_raysamples = 1 + np.ceil(MIN_SAMPLE_FACTOR / DX * np.abs(z) * np.max((np.abs((np.tan(theta1)-np.tan(theta2))), np.abs((np.tan(phi1)-np.tan(phi2)))) )); # this was set to "10" before
                        if (u == 0 and v==0) and False: # or (u == self.nu/2 and v == self.nv/2):
                            print 'z,u,v,num_raysamples, num_zsamples',z,u,v,num_raysamples, num_zsamples
                    num_zsamples = np.int(num_zsamples)
                    num_raysamples = np.int(num_raysamples)
                    from pylflib import compute_rayspread
                    spread = compute_rayspread(self, u, v, 0, z_spacing, dummy_radiometry,
                                               sine_correction = True,
                                               supersample_factor = supersample_factor, num_raysamples = num_raysamples,
                                               num_zsamples =num_zsamples, pixel_fill_factor = 0.9, ulens_fill_factor =ulens_fill_factor)
                    psf = spread[0]

                    # compute base radiometry
                    radiometry[v,u] = psf.sum()
        
        # ---------------------
        # Compute Ray Spreads
        # ---------------------
        import time
        tic = time.time()

        # Compute the rayspreads
        rayspread_db = RayspreadDatabase(self.nu, self.nv, self.ns, self.nt, z_center, z_coords, supersample_factor)
        coords = []
        for u in range(self.nu):
            for v in range(self.nv):
                for z_index in range(z_len):
                    z = z_coords[z_index]
                    
                    
                    #optimize = False
                    if not optimize:
                        # default sampling values                    
                        num_raysamples = 10
                        num_zsamples = 5
                    else:
                        # optimize raypspread computation by having a z-depth and angle dependent sampling rate
                        # (see Evernote 'rayspread code optimization')
                        #(theta, phi) = self.uv_to_angle(u, v)
                        DX = self.ulens_pitch * 1.0 / self.objective_magnification / supersample_factor;
                        MIN_ZSAMPLE_FACTOR = 1.0; # increase this for finer z-sampling
                        THETA_NA = np.arcsin(self.objective_na / self.medium_index);
                        num_zsamples = np.ceil(z_spacing * MIN_ZSAMPLE_FACTOR / DX * np.tan(THETA_NA)); # this was set to "5" before
                        
                        (theta1, phi1) =   self.uv_to_angle( u-pixel_fill_factor/2, v-pixel_fill_factor/2) # this is only approximate, since this uses the old incorrect uv_to_angle()
                        (theta2, phi2) =   self.uv_to_angle( u+pixel_fill_factor/2, v+pixel_fill_factor/2)
                        MIN_SAMPLE_FACTOR = 1.0; # increase this for finer (u,v) sampling
                        num_raysamples = 1 + np.ceil(MIN_SAMPLE_FACTOR / DX * np.abs(z) * np.max((np.abs((np.tan(theta1)-np.tan(theta2))), np.abs((np.tan(phi1)-np.tan(phi2)))) )); # this was set to "10" before
                        if (u == 0 and v==0) and False: # or (u == self.nu/2 and v == self.nv/2):
                            print 'z,u,v,num_raysamples, num_zsamples',z,u,v,num_raysamples, num_zsamples
                    num_zsamples = np.int(num_zsamples)
                    num_raysamples = np.int(num_raysamples)
                    coords.append((self, u, v, z, z_spacing, radiometry, supersample_factor, z_index,num_raysamples, num_zsamples, pixel_fill_factor, ulens_fill_factor, ulens_profile))
        if optimize and False:
            print 'DX, THETA_NA', DX, THETA_NA * 180 / np.pi
        print  'len(coords)',len(coords)
        #raise SystemExit #TEMP
        print '\tStarting multiprocessing pool with ', num_threads , ' worker processes.'
        from multiprocessing import Pool
        pool = Pool(processes=num_threads)
        results = pool.imap(compute_rayspread_imap, coords, chunksize=100)
        for r in results:

            # Extract the data as it streams in from the worker
            # threads, and print out a progress message.
            out_tuple = r[0]
            if out_tuple[2] == 0 and out_tuple[7] == 0:
                print '\tGenerating rayspreads. (%0.2f%% complete)' % (float(out_tuple[1]+1)/rayspread_db.nu*100)

            u = out_tuple[1]
            v = out_tuple[2]
            z_index = out_tuple[7]


            spread = r[1]
            psf = spread[0]
            kernel_width = psf.shape[1]
            kernel_height = psf.shape[0]
            kernel_h_offset = spread[1]
            kernel_v_offset = spread[2]

            # Degenerate rayspreads will have zero width and height.  We ignore these.
            if kernel_width == 0 or kernel_height == 0:
                continue

            # This is a useful debugging check that helps us to avoid buggy rayspreads.
            if psf.min() != psf.min():
                print 'ERROR: The PSF contains a NaN value for rayspread (%d, %d @ z=%d).  These rayspreads need further debugging!' % (u, v, z_index)
                exit(1)

            # Generate the flipped PSF (deprecated: only used by the CPU project() and backproject() methods)
            flipped_psf = np.fliplr(np.flipud(spread[0]))
            flipped_kernel_width = kernel_width
            flipped_kernel_height = kernel_height
            flipped_kernel_h_offset = -(kernel_h_offset + kernel_width - 1)
            flipped_kernel_v_offset = -(kernel_v_offset + kernel_height - 1)

            # Store the result
            rayspread_db.rayspreads[(z_index,u,v)] = (psf, kernel_h_offset, kernel_v_offset)
            rayspread_db.flipped_rayspreads[(z_index,u,v)] = (flipped_psf,
                                                              flipped_kernel_h_offset,
                                                              flipped_kernel_v_offset)
            

        pool.close()
        pool.join()
        print 'rayspread calculation took', time.time() - tic ,'seconds.'

        return rayspread_db        

    # --------------------
    # WAVESPREAD GENERATION
    # --------------------
    def trim_to_na(self, coords, coefficients):
        '''
        Given a set of coefficients in a light field splat function, this
        routine removes any coefficients that are outside the numerical
        aperture of the objective.
        '''

        # DEBUGGING FLAGS
        DEBUG = 0

        # Step 1: Compute padding.  This depends on optical parameters 
        alpha_o = np.arcsin(self.objective_na / self.medium_index)
        objective_focal_length = self.tubelens_focal_length / self.objective_magnification
        x_o = objective_focal_length * self.medium_index * np.sin(alpha_o)

        alpha_i = np.arctan2( self.ulens_pitch/2.0, self.ulens_focal_length )
        x_i = self.tubelens_focal_length * np.tan(alpha_i)

        trim_padding = int(np.ceil(self.nu * ( 1.0 - x_o/x_i )))

        if DEBUG:
            print '--> Trimming to objective\'s numerical aperture'
            print '    alpha_object           : %0.2f' % (180.0 / np.pi * alpha_o)
            print '    alpha_image            : %0.2f' % (180.0 / np.pi * alpha_i)
            print '    back aperture diameter : %0.2f' % (x_o)
            print '    lenslet ap. diameter   : %0.2f' % (x_i)
            print '    trim_padding            : %d pixels' % (trim_padding)

        u_coords = coords[:,0]
        v_coords = coords[:,1]
        entries_within_na = np.nonzero(np.sqrt(np.power(u_coords-(self.nu-1.0)/2.0,2) +
                                               np.power(v_coords-(self.nv-1.0)/2.0,2)) < (self.nv/2.0 - trim_padding))[0]
        before_count = coords.shape[0]
        coefficients = coefficients[entries_within_na]
        coords = coords[entries_within_na,:]
        if DEBUG:
            print '    trimmed %d / %d coefficients ( %0.2f%% )' % (before_count,
                                                                    coords.shape[0],
                                                                    float(coords.shape[0]) / before_count * 100)

        return (coords, coefficients)

    
    def physical_rayspread_database(self, z_len, z_center, z_spacing, radiometry=None,
                                    supersample_factor = 1, voxels_as_points = False,
                                    num_threads = 1, optimize = True):

        """
        Generate a rayspread database for this lenslet array based on
        the physical optics model.  This routine uses the Arroyo C++
        library to produce a light field PSF that accounts for
        diffraction.  This database in needed by the reconstruction
        and projection algorithms in volume.py.

        z_len             - the number of slices in the stack
        z_center          - the z offset in micron from the focal plane
                            of the center slice in the volume
                            positive z is "up" from the focal plane towards
                            the objective.  If there are an even number of
                            slices in the volume, then this is the average
                            of the z offset of the two middle slices.
        z_spacing         - the number of micron between slices in the
                            volume
        radiometry        - is either 'None' or a mapping (u,v)->float 
                            that specifies the intensity of a uniform
                            fluorescent volume
                            along each of the ray angles.  All the (u,v)
                            ray angles that will be used in the reconstruction
                            must be specified in this mapping.
        supersample_factor - Request supersampled ray spread kernels.  This is
                             useful in super-resolved deconvolution.
        voxels_as_points   - If True, treat each voxel as an ideal point source.
                             This turns of numerical integration that gives the
                             voxel spatial extent (which can be important for anti-aliasing.
        """
        import time
        tic = time.time()

        # modify these to adjust the fill factor used by the rayspread code
        pixel_fill_factor = self.pixel_fill_factor
        ulens_fill_factor = self.ulens_fill_factor
        effective_ulens_pitch = self.ulens_pitch / self.objective_magnification

        # ---------------------
        # Compute z coordinates
        # ---------------------
        #
        # The coordinates are listed from the bottom slice in the
        # stack (slice 0) to the top.  The convention is that +z is
        # closer to the objective lens in object space.
        #
        z_start = z_center - 0.5 * (z_len - 1) * z_spacing
        z_coords = [i * z_spacing + z_start for i in range(z_len)]

        # ----------------
        # Generate LF psfs 
        # ----------------
        # 
        # x and y coords are voxel centers covering a single lenslet
        # When supersample_factor is 1, this just yields (0,0)
        x_coords = np.linspace(-0.5+1.0/(2*supersample_factor),
                                0.5-1.0/(2*supersample_factor),
                                supersample_factor) * effective_ulens_pitch
        y_coords = np.linspace(-0.5+1.0/(2*supersample_factor),
                                0.5-1.0/(2*supersample_factor),
                                supersample_factor) * effective_ulens_pitch
        print '\tX Coords:', x_coords
        print '\tY Coords:', y_coords
        print '\tZ Coords:', z_coords

        # The ASPF cache stores previously computing APSFs for specific
        # z-depths.  When building super-sampled wavespreads, these APSFs are
        # re-used several times so it saves considerable time if we store and
        # re-use the APSF lookup tables.
        apsf_cache = ApsfCache()

        # Compute the point spread functions
        psf_db = PsfDatabase(self.nu, self.nv, self.ns, self.nt, z_center, z_coords, supersample_factor)
        psf_db.x_coords = x_coords
        psf_db.y_coords = y_coords
        psf_db.z_coords = z_coords
        
        imap_args = []
        for x in range(len(x_coords)):
            for y in range(len(y_coords)):
                for z in range(len(z_coords)):
                    objective_theta = np.arcsin(self.objective_na / self.medium_index);
                    aperture_diameter = np.abs( 2*z_coords[z]*np.tan(objective_theta) );

                    num_lenslets_in_aperture = int(np.ceil(aperture_diameter /
                                                           (self.ulens_pitch / self.objective_magnification)))

                    # Make sure that there are an odd number of lenslets in the aperture.
                    if num_lenslets_in_aperture % 2 == 0:
                        num_lenslets_in_aperture += 1

                    # Add some extra padding...
                    num_lenslets_in_aperture += 2

                    imap_args.append( (x_coords[x] * 1e-6, y_coords[y] * 1e-6, z_coords[z] * 1e-6,
                                       z_spacing, supersample_factor, 
                                       num_lenslets_in_aperture, self, radiometry, x, y, z, voxels_as_points,
                                       apsf_cache) )

                    # Debugging
                    # compute_psf_imap(imap_args[0])
                    # raise SystemExit

        print '\tStarting multiprocessing pool with ', num_threads , ' worker processes.'
        from multiprocessing import Pool
        pool = Pool(processes=num_threads)
        results = pool.imap(compute_psf_imap, imap_args, chunksize=1)
        for r in results:

            x = r[2]
            y = r[3]
            z = r[4]
            num_lenslets_in_aperture = r[5]
            psf_sum = r[6]
            num_nonzero = r[7]
            pct_pixels_retained = r[8]
            pct_intensity_retained = r[9]


            # Store the coordinates and coefficients as
            # matrices in the psf_db data structure.  Each row
            # of these matrices contains one non-zero data
            # element from this psf.
            psf_db.psf_coordinates[(x,y,z)] = r[0]
            psf_db.psf_coefficients[(x,y,z)] = r[1]
            psf_db.max_coefficient = np.maximum(psf_db.max_coefficient, r[1].max())

            # Print status update
            print '\t[%d %d] @ z = %d (%f um) -- psf diameter = %d lenslets  sum = %0.2f.  Kept %d / %d  nonzero coefficients (%0.1f%% pixels,  %0.2f%% intensity)' % (x, y, z, z_coords[z], num_lenslets_in_aperture, psf_sum, r[1].shape[0], num_nonzero, pct_pixels_retained, pct_intensity_retained)


            
                                        
        print 'light field psf calculation took', time.time() - tic ,'seconds.'
        return psf_db
