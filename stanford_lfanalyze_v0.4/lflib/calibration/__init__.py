import math, datetime
import numpy as np
from lflib.lfexceptions import ZeroImageException
from lflib.lightfield import LightField
from lflib.calibration.imaging import CalibrationAlignmentMethods

class LightFieldCalibration(object):


    def __init__(self, ulens_focal_length, ulens_focal_distance, ulens_pitch, pixel_size,
                 objective_magnification, objective_na, medium_index,
                 tubelens_focal_length,
                 ulens_fill_factor, pixel_fill_factor, ulens_profile, center_wavelength,
                 alignment_method = CalibrationAlignmentMethods.CALIBRATION_ISOMETRY_ALIGNMENT):

        # Optical parameters
        self.ulens_focal_length = ulens_focal_length
        self.ulens_focal_distance = ulens_focal_distance
        self.ulens_pitch = ulens_pitch
        self.pixel_size = pixel_size
        self.objective_magnification = objective_magnification
        self.objective_na = objective_na
        self.medium_index = medium_index
        self.tubelens_focal_length = tubelens_focal_length
        self.ulens_fill_factor = ulens_fill_factor
        self.pixel_fill_factor = pixel_fill_factor
        self.ulens_profile = ulens_profile
        self.center_wavelength = center_wavelength
        
        # Geometric Calibration
        self.geometry_is_calibrated = False
        self.calibration_image_filename = None
        self.skip_alignment = False
        self.alignment_method = alignment_method
        self.geometric_calibration_model = None
        self.radiometric_calibration_model = None
        self.nu = None
        self.nv = None
        self.ns = None
        self.nt = None
        
        # Radiometric Calibration
        self.radiometry_is_calibrated = False
        self.has_dark_frame = False
        self.dark_frame_filename = None
        self.dark_frame = None
        self.radiometric_correction = None
        self.rectified_radiometry = None
        self.aperture_radiometry = None
        self.radiometry = None
        
        # Rayspread database
        self.rayspread_db_generated = False
        self.rayspread_db = None
        self.psf_db = None
        self.lenslet_array = None
        
        # - Zahid: Adding volume calibration too
        self.num_z_slices   = None
        self.um_per_z_slice = None
        self.z_center     = None

    def compute_trim_padding(self, printout = False):
        alpha_o = np.arcsin(self.objective_na / self.medium_index)
        objective_focal_length = self.tubelens_focal_length / self.objective_magnification
        x_o = objective_focal_length * self.medium_index * np.sin(alpha_o)

        alpha_i = np.arctan2( self.ulens_pitch/2.0, self.ulens_focal_length )
        x_i = self.tubelens_focal_length * np.tan(alpha_i)

        trim_padding = int(np.ceil(self.nu * ( 1.0 - x_o/x_i )))

        if printout:
            print 'Objective / lenslet aperture information:'
            print '    alpha_object           : %0.2f degrees' % (180.0 / np.pi * alpha_o)
            print '    alpha_image            : %0.2f degrees' % (180.0 / np.pi * alpha_i)
            print '    back aperture diameter : %0.2f mm' % (x_o)
            print '    lenslet ap. diameter   : %0.2f mm' % (x_i)
            print '    trim_padding           : %d pixels' % (trim_padding)

        return trim_padding


    def mask_trimmed_angles(self, lightfield):
        lf_subaperture = lightfield.asimage(representation = LightField.TILED_SUBAPERTURE)

        trim_padding = self.compute_trim_padding()
        nCols = lf_subaperture.shape[1] / self.nu
        nRows = lf_subaperture.shape[0] / self.nv

        for u in range(self.nu):
            for v in range(self.nv):
                r = np.sqrt(np.power(u-(self.nu-1.0)/2.0,2) + np.power(v-(self.nv-1.0)/2.0,2))
                if r >= (self.nv/2.0 - trim_padding):
                    lf_subaperture[v*nRows:v*nRows+nRows, u*nCols:u*nCols+nCols] = 0

        return LightField(lf_subaperture, self.nu, self.nv, self.ns, self.nt,
                          representation = LightField.TILED_SUBAPERTURE)
        
    
    def print_summary(self):
        print '\n------------------ LF Calibration Report --------------------'
        print 'ulens_focal_length (um)  : ', self.ulens_focal_length
        print 'ulens_focal_distance (um): ', self.ulens_focal_distance
        print 'ulens_pitch  (um)        : ', self.ulens_pitch
        print 'ulens_profile            : ', self.ulens_profile
        print 'pixel_size   (um)        : ', self.pixel_size
        print 'objective_magnification  : ', self.objective_magnification
        print 'objective_na             : ', self.objective_na
        print 'medium_index             : ', self.medium_index
        print 'tubelens_focal_length    : ', self.tubelens_focal_length

        # Print NA trimming info
        self.compute_trim_padding(printout = True)
        
        print 'ulens_fill_factor        : ', self.ulens_fill_factor, '\n'
        print 'pixel_fill_factor        : ', self.pixel_fill_factor, '\n'
        print '(nu,nv,ns,nt)            : ', self.nu, self.nv, self.ns, self.nt
        print '\nGeometry calibrated    : ', self.geometry_is_calibrated
        if self.geometry_is_calibrated and self.geometric_calibration_model != None:
            print 'Geometric calibration'
            print 'forward : \n', self.geometric_calibration_model.forwardCoefficients
            print 'reverse : \n', self.geometric_calibration_model.reverseCoefficients
        if self.geometry_is_calibrated and self.radiometric_calibration_model != None:
            print 'Geometric calibration for radiometry image'
            print '\tforward : \n', self.radiometric_calibration_model.forwardCoefficients
            print '\treverse : \n', self.radiometric_calibration_model.reverseCoefficients
        print '\nRadiometry calibrated  : ', self.radiometry_is_calibrated
        if self.radiometry_is_calibrated:
            print '\tMin %f   Max: %f   Mean: %f' % (self.radiometry.min(), self.radiometry.max(), self.radiometry.mean())
            print '\tDarkframe:', self.has_dark_frame
        print '\nRayspread DB generated : ', self.rayspread_db_generated
        if self.rayspread_db_generated:
            print '\tnumber of z slices   : ',self.num_z_slices
            print '\tspacing per z slices : ',self.um_per_z_slice
            print '\tz center             : ',self.z_center
            if self.rayspread_db:
                print '\trayspreads in the database: %d' % (len(self.rayspread_db.rayspreads.keys()))
            if self.psf_db:
                print '\tlight field psfs in the database: %d' % (len(self.psf_db.psf_coordinates.keys()))
        print '-------------------------------------------------------------'
    
    # ----------------------- SAVE/LOAD ROUTINES -----------------------------
    
    @classmethod
    def load(cls, filename):
        '''
        This method loads calibration data from a file on disk.  You can use it like this:
          
          calibration_instance = LightFieldCalibration.load(filename)
        '''
        import h5py
        calibration_file = h5py.File(filename, 'r')
        
        # Optics
        ulens_focal_length = calibration_file['/optics'].attrs['ulens_focal_length']
        try:  # backwards compatibility with older calibration files
            ulens_focal_distance = calibration_file['/optics'].attrs['ulens_focal_distance']
        except KeyError, e:
            ulens_focal_distance = ulens_focal_length
        ulens_pitch = calibration_file['/optics'].attrs['ulens_pitch']
        pixel_size = calibration_file['/optics'].attrs['pixel_size']
        objective_na = calibration_file['/optics'].attrs['objective_na']
        objective_magnification = calibration_file['/optics'].attrs['objective_magnification']
        medium_index = calibration_file['/optics'].attrs['medium_index']
        tubelens_focal_length = calibration_file['/optics'].attrs['tubelens_focal_length']
        ulens_fill_factor = calibration_file['/optics'].attrs['ulens_fill_factor']
        pixel_fill_factor = calibration_file['/optics'].attrs['pixel_fill_factor']
        center_wavelength = calibration_file['/optics'].attrs['center_wavelength']
        ulens_profile = calibration_file['/optics'].attrs['ulens_profile']
        
        # Create the class instance.
        instance = cls(ulens_focal_length, ulens_focal_distance, ulens_pitch, pixel_size,
                       objective_magnification, objective_na, medium_index,
                       tubelens_focal_length,
                       ulens_fill_factor, pixel_fill_factor,
                       ulens_profile, center_wavelength)
        
        # Geometry
        try:
            geometry_group = calibration_file['/geometry']
            instance.calibration_image_filename = geometry_group.attrs['calibration_image_filename']
            instance.skip_alignment = geometry_group.attrs['skip_alignment']
            instance.alignment_method = geometry_group.attrs['alignment_method']
            instance.nu = geometry_group.attrs['nu']
            instance.nv = geometry_group.attrs['nv']
            instance.ns = geometry_group.attrs['ns']
            instance.nt = geometry_group.attrs['nt']
            
            # TODO: It would be nice if image warp models serialized
            # themselves somehow so that we could stor them in a more
            # generic manner.  For now we directly store affine
            # coefficients.
            if not instance.skip_alignment:
                from lflib.calibration import models
                from lflib.calibration.imaging import CalibrationAlignmentMethods
                if instance.alignment_method == CalibrationAlignmentMethods.CALIBRATION_AFFINE_ALIGNMENT:
                    instance.geometric_calibration_model = models.AffineWarp()
                    instance.radiometric_calibration_model = models.AffineWarp()
                elif instance.alignment_method == CalibrationAlignmentMethods.CALIBRATION_CUBIC_ALIGNMENT:
                    instance.geometric_calibration_model = models.CubicWarp()
                    instance.radiometric_calibration_model = models.CubicWarp()
                else:
                    instance.geometric_calibration_model = models.IsometryWarp()
                    instance.radiometric_calibration_model = models.IsometryWarp()
                    
                instance.geometric_calibration_model.forwardCoefficients = geometry_group['forwardCoefficients'][...]
                instance.geometric_calibration_model.reverseCoefficients = geometry_group['reverseCoefficients'][...]

                instance.radiometric_calibration_model.forwardCoefficients = geometry_group['radiometryForwardCoefficients'][...]
                instance.radiometric_calibration_model.reverseCoefficients = geometry_group['radiometryReverseCoefficients'][...]
            
            instance.geometry_is_calibrated = True
        except KeyError, e:
            instance.geometry_is_calibrated = False
        
        # Radiometry
        try:
            radiometry_group = calibration_file['radiometry']
            instance.radiometry = radiometry_group['radiometry'][...]
            instance.radiometric_correction = radiometry_group['radiometric_correction'][...]
            instance.rectified_radiometry = radiometry_group['rectified_radiometry'][...]
            instance.has_dark_frame = radiometry_group.attrs['has_dark_frame']
            if instance.has_dark_frame:
                instance.dark_frame_filename = radiometry_group.attrs['dark_frame_filename']
                instance.dark_frame = radiometry_group['dark_frame'][...]
            
            instance.radiometry_is_calibrated = True
        except KeyError, e:
            instance.has_dark_frame = False
            instance.radiometry_is_calibrated = False
        
        # Ray spreads
        try:
            rayspread_group = calibration_file['rayspread']
            psf_group = calibration_file['psf']
            
            instance.num_z_slices   = rayspread_group.attrs['num_z_slices']
            instance.um_per_z_slice = rayspread_group.attrs['um_per_z_slice']
            instance.z_center       = rayspread_group.attrs['z_center']

            from lflib.optics import LensletArray
            instance.lenslet_array = LensletArray(instance.nu, instance.nv, instance.ns, instance.nt,
                                                  instance.ulens_pitch, instance.pixel_size,
                                                  instance.ulens_focal_length, instance.ulens_focal_distance,
                                                  instance.objective_magnification, instance.objective_na,
                                                  instance.medium_index, instance.tubelens_focal_length,
                                                  instance.ulens_fill_factor, instance.pixel_fill_factor,
                                                  instance.ulens_profile, instance.center_wavelength)

            import pickle
            try:
                instance.rayspread_db = pickle.loads(rayspread_group['rayspread_db'][...].tostring())
            except:
                import traceback
                #print 'Caught exception loading rayspread_db:' 
                #print traceback.format_exc()
                instance.rayspread_db = None
            try:
                instance.psf_db = pickle.loads(psf_group['psf_db'][...].tostring())
            except:
                import traceback
                #print 'Caught exception loading psf_db:' 
                #print traceback.format_exc()
                instance.psf_db = None
                
            instance.rayspread_db_generated = True

        except KeyError, e:
            print e
            instance.rayspread_db_generated = False
        
        calibration_file.close()
        return instance

    
    def save(self, filename):
        '''
        Save the calibartion image to disk.  Calibration information
        is stored in an hdf5 file using attributes.
        '''
        import h5py
        calibration_file = h5py.File(filename, 'w')
        calibration_file.attrs['timestamp'] = str(datetime.datetime.now())
        calibration_file.attrs['comments'] = ''    # Empty for now...
        
        # Optics
        optics_group = calibration_file.create_group('optics')
        optics_group.attrs['ulens_focal_length'] = self.ulens_focal_length
        optics_group.attrs['ulens_focal_distance'] = self.ulens_focal_distance
        optics_group.attrs['ulens_pitch'] = self.ulens_pitch
        optics_group.attrs['pixel_size'] = self.pixel_size
        optics_group.attrs['objective_magnification'] = self.objective_magnification
        optics_group.attrs['objective_na'] = self.objective_na
        optics_group.attrs['medium_index'] = self.medium_index
        optics_group.attrs['tubelens_focal_length'] = self.tubelens_focal_length
        optics_group.attrs['ulens_fill_factor'] = self.ulens_fill_factor
        optics_group.attrs['pixel_fill_factor'] = self.pixel_fill_factor
        optics_group.attrs['center_wavelength'] = self.center_wavelength
        optics_group.attrs['ulens_profile'] = self.ulens_profile
        
        # Geometric Calibration
        if self.geometry_is_calibrated:
            geometry_group = calibration_file.create_group('geometry')
            geometry_group.attrs['calibration_image_filename'] = self.calibration_image_filename
            geometry_group.attrs['skip_alignment'] = self.skip_alignment
            geometry_group.attrs['alignment_method'] = self.alignment_method
            geometry_group.attrs['nu'] = self.nu
            geometry_group.attrs['nv'] = self.nv
            geometry_group.attrs['ns'] = self.ns
            geometry_group.attrs['nt'] = self.nt
            
            # TODO: It would be nice if image warp models serialized
            # themselves somehow so that we could stor them in a more
            # generic manner.  For now we directly store the
            # coefficient matrix.
            if not self.skip_alignment:
                geometry_group.create_dataset('forwardCoefficients', data = self.geometric_calibration_model.forwardCoefficients)
                geometry_group.create_dataset('reverseCoefficients', data = self.geometric_calibration_model.reverseCoefficients)

                # If a radiometry image was supplied, it will have its
                # own geometry calibration.  Otherwise is shares the
                # default geometric calibration.
                if self.radiometric_calibration_model != None:
                    geometry_group.create_dataset('radiometryForwardCoefficients',
                                                  data = self.radiometric_calibration_model.forwardCoefficients)
                    geometry_group.create_dataset('radiometryReverseCoefficients',
                                                  data = self.radiometric_calibration_model.reverseCoefficients)
                else:
                    geometry_group.create_dataset('radiometryForwardCoefficients',
                                                  data = self.geometric_calibration_model.forwardCoefficients)
                    geometry_group.create_dataset('radiometryReverseCoefficients',
                                                  data = self.geometric_calibration_model.reverseCoefficients)
            
        if self.radiometry_is_calibrated:
            radiometry_group = calibration_file.create_group('radiometry')
            radiometry_group.create_dataset('radiometry', data = self.radiometry)
            radiometry_group.create_dataset('radiometric_correction', data = self.radiometric_correction)
            radiometry_group.create_dataset('rectified_radiometry', data = self.rectified_radiometry)
            radiometry_group.attrs['has_dark_frame'] = self.has_dark_frame
            if self.has_dark_frame:
                radiometry_group.attrs['dark_frame_filename'] = self.dark_frame_filename
                radiometry_group.create_dataset('dark_frame', data = self.dark_frame, compression = 'gzip')

        
        # Rayspreads
        if self.rayspread_db_generated:
            print 'Rayspread db generated'
            rayspread_group = calibration_file.create_group('rayspread')
            psf_group = calibration_file.create_group('psf')
            
            # saving out the volume calibration data too
            rayspread_group.attrs['num_z_slices']   = self.num_z_slices
            rayspread_group.attrs['um_per_z_slice'] = self.um_per_z_slice
            rayspread_group.attrs['z_center']       = self.z_center

            # TODO: Pickling the rayspread_db is definitely the
            # fastest way to save and load it for now, but it's
            # brittle and opaque, and won't be backwards compatible if
            # we ever change the rayspread_db class.
            import pickle
            if self.rayspread_db != None:
                rayspread_group.create_dataset('rayspread_db',
                                               data = np.fromstring(pickle.dumps(self.rayspread_db),
                                                                    dtype='uint8'),
                                               compression = 'gzip')

            if self.psf_db != None:
                psf_group.create_dataset('psf_db',
                                         data = np.fromstring(pickle.dumps(self.psf_db), dtype='uint8'),
                                         compression = 'gzip')
        
        calibration_file.close()

    
    # ----------------------- GEOMETRY --------------------------------
    
    def calibrate_geometry(self, calibration_image_file, skip_alignment = False,
                           skip_subpixel_alignment = False,
                           debug_mode = False, chief_ray_image = False,
                           radiometry_file = None, align_radiometry = False,
                           crop_center_lenslets = False, lenslet_detect_thres=0.3):
        '''
        Run the light field imaging geometric calibration routine.
        '''
        
        from lflib.imageio import load_image
        self.calibration_image_filename = calibration_image_file
        calibration_image = load_image(calibration_image_file, dtype=np.float32, normalize = True)

        # If the image is RGB, we average the colors here.  This isn't the
        # best grayscale conversion, but it certainly will suffice for
        # calibration purposes.
        if len(calibration_image.shape) == 3:
            calibration_image = np.mean(calibration_image, axis=2)
        
        # Do a quick check to see if we capture an all-black (or mostly
        # black) frame.  If so, we should probably alert the user!!
        if (calibration_image.mean() < 0.005 and not skip_alignment ):
            raise Exception("Error: calibration image was too dark, mean=%f."%calibration_image.mean() +
                            " Please check input settings and try again.")
        
        # Step 1: run the calibration routine
        self.skip_alignment = skip_alignment
        if self.skip_alignment:
            self.geometric_calibration_model = None
            approx_ppl = int(np.ceil(float(self.ulens_pitch) / self.pixel_size))
            self.nu = approx_ppl
            self.nv = approx_ppl
            self.ns = calibration_image.shape[1] / self.nu
            self.nt = calibration_image.shape[0] / self.nv
        else:
            print 'Performing geometric calibration with: ', calibration_image_file

            from lflib.calibration.imaging import geometric_calibration
            try: 
                (self.geometric_calibration_model,
                 self.geometric_calibration_residuals) = geometric_calibration(calibration_image,
                                                                               self.pixel_size, self.ulens_pitch,
                                                                               self.alignment_method,
                                                                               chief_ray_image,
                                                                               crop_center_lenslets)
            except RuntimeError, e:
                print str(e)
                print '\nPlease check your calibration parameters and try again.'
                raise SystemExit
            
            if radiometry_file != None and align_radiometry == True:
                print 'Performing a separate geometric calibration on radiometry frame:', radiometry_file
                radiometry_image = load_image(radiometry_file, dtype=np.float32, normalize = True)

                # If the image is RGB, we average the colors here.  This isn't the
                # best grayscale conversion, but it certainly will suffice for
                # calibration purposes.
                if len(radiometry_image.shape) == 3:
                    radiometry_image = np.mean(radiometry_image, axis=2)

                # Do a quick check to see if we capture an all-black (or mostly
                # black) frame.  If so, we should probably alert the user!!
                if ( radiometry_image.mean() < 0.01 and not self.skip_alignment ):
                    raise Exception("Error: radiometry image was too dark. " +
                                    " Please check input settings and try again.")

                try:
                    (self.radiometric_calibration_model,
                     self.radiometric_calibration_residuals) = geometric_calibration(radiometry_image,
                                                                                     self.pixel_size,
                                                                                     self.ulens_pitch, 
                                                                                     self.alignment_method,
                                                                          False, # No chief ray calibration
                                                                                     crop_center_lenslets)
                except RuntimeError, e:
                    print str(e)
                    print '\nPlease check your calibration parameters and try again.'
                    raise SystemExit

            else:
                self.radiometric_calibration_model = self.geometric_calibration_model
                self.radiometric_calibration_residuals = self.geometric_calibration_residuals

            # Step 2: compute the number of spatial and ray angle 
            # samples in the calibrated light field (ns, nt, nu, nv).
            output_pixels_per_lenslet = int(math.ceil( float(self.ulens_pitch) / self.pixel_size ))
            self.nu = output_pixels_per_lenslet
            self.nv = output_pixels_per_lenslet

            rectified_im = self.geometric_calibration_model.warp_image(calibration_image,
                                                                       output_pixels_per_lenslet,
                                                                       'r', cropToInside = True)
            
            self.ns = rectified_im.shape[1]/self.nu
            self.nt = rectified_im.shape[0]/self.nv
        
        # Set the flag that geometry is now calibrated
        self.geometry_is_calibrated = True
    
    def rectify_lf(self, light_field_image, model = None, cropToInside = True,
                   lenslet_offset = None, output_size = None):
        '''
        Rectify the image.  The routine automatically selects a
        suitable output size so that the entire image contains valid,
        whole lenslets.
        '''
        if model == None:
            model = self.geometric_calibration_model

        rectified_im = model.warp_image(light_field_image,
                                        max(self.nu, self.nv), 'r',
                                        cropToInside = cropToInside,
                                        lenslet_offset = lenslet_offset,
                                        output_size = output_size)

        return LightField(rectified_im, self.nu, self.nv, self.ns, self.nt,
                          representation = LightField.TILED_LENSLET)

    
    # ----------------------- RADIOMETRY --------------------------------
    def set_radiometry(self, radiometry, mask_to_na = True):
        self.radiometry = radiometry
        
        # Drop the pixels that are outside of the numerical aperture of
        # the objective and the lenslet array.  These "corner" rays do
        # little help with the reconstruction; they mostly end up adding
        # stray light and noise.
        if mask_to_na:
            nv = self.radiometry.shape[0]
            nu = self.radiometry.shape[1]
            for u in range(nu):
                for v in range(nv):
                    if (np.sqrt(np.power(u-(nu-1.0)/2.0,2) + np.power(v-(nv-1.0)/2.0,2)) >= nv/2.0):
                        self.radiometry[v,u] = 0.0;
    
    def get_radiometry(self):
        return self.radiometry

    def correct_radiometry(self, im):
        if self.radiometric_correction != None:
            print '\t--> Applying radiometric correction using the full model'
            result = im * self.radiometric_correction
        else:
            print '\tWARNING: no radiometric calibration found.  Skipping.'
            result = im  # If no radiometric correction exists in this calibration object, skip it!
        return LightField(result, self.nu, self.nv, self.ns, self.nt,
                          representation = LightField.TILED_LENSLET)


    def subtract_dark_frame(self, raw_image):
        '''
        Subtract the dark frame from a raw (unrectified) input image.
        If no dark frame exists for this calibration, the same image
        is returned.
        '''
        if self.has_dark_frame:
            input_dtype = raw_image.dtype
            output_im = raw_image.astype(np.float32) - self.dark_frame

            if output_im.max() <= 0.0:
                raise ZeroImageException("ERROR: dark frame subtracted image" 
                        "contains only zeros. Perhaps you are using the wrong" 
                        "dark frame?")
            
            return output_im.astype(input_dtype)
        else:
            return raw_image
    
    def calibrate_radiometry(self, calibration_image_file, radiometry_frame_file = None, dark_frame_file = None):
        
        from lflib.imageio import load_image, save_image
        
        if radiometry_frame_file != None:
            im = load_image(radiometry_frame_file, dtype=np.float32, normalize = False)
        else:
            im = load_image(calibration_image_file, dtype=np.float32, normalize = False)
        
        if dark_frame_file is not None:
            print '\t--> Subtracting dark frame: ', dark_frame_file
            self.dark_frame_filename = dark_frame_file
            self.dark_frame = load_image(dark_frame_file, dtype=np.float32, normalize = False)
            self.has_dark_frame = True
            im -= self.dark_frame

            if im.max() <= 0.0:
                raise ZeroImageException("ERROR: image contains only zeros.")

        if self.radiometric_calibration_model != None:
            print '\t--> rectifying using radiometry-specific geometric calibration.'

            # The calibration routine may align the geometry and
            # radiometry images to each other with some integer number
            # lenslet shift between them.  We assume that the
            # radiometry image should be within one lenslet of the
            # geometry image, so we compensate for this here.
            radcal_upper_left = self.radiometric_calibration_model.eval_point([0.0,0.0], 'f')
            geomcal_upper_left = self.geometric_calibration_model.eval_point([0.5,0.5], 'f')
            diff = (radcal_upper_left - geomcal_upper_left) / 0.5
            shift = np.multiply(np.floor(np.abs(diff)), np.sign(diff))
            if shift.sum() != 0.0: 
                print '\t    Detected a lenslet shift between the geometry and radiometry images:', shift, ' Compensating.'

            rectified_lf = self.rectify_lf(im, self.radiometric_calibration_model,
                                           lenslet_offset = shift)

            # Disabling this for now... may still be required when --align-radiometry is used... -mbroxton
                                           #output_size = (self.nv*self.nt, self.nu*self.ns))
        elif not self.skip_alignment:
            print '\t--> rectifying using default geometric calibration.'
            rectified_lf = self.rectify_lf(im, self.geometric_calibration_model)
        else:
            rectified_lf = LightField(im, self.nu, self.nv, self.ns, self.nt,
                                      representation = LightField.TILED_LENSLET)
        from lflib.imageio import save_image

        from lflib.volume import LightFieldProjection
        lfproj = LightFieldProjection(self.rayspread_db, self.psf_db, gpu_id = 1)
        if self.psf_db:
            vol_ones = self.psf_db.ones_volume()
        else:
            vol_ones = self.rayspread_db.ones_volume()
        self.ideal_lf = lfproj.project(vol_ones).asimage(representation = LightField.TILED_LENSLET)

        # "BASIC" RADIOMETRY: Compute the average aperture vignetting function
        #
        # If no radiometry image is provided, we compute an
        # approximate vignetting function by averaging together all of
        # the data for each direction in the aperture.
        aperture_avg = np.zeros([self.nv, self.nu], dtype=np.float32);
        for v in range(self.nv):
            for u in range(self.nu):
                subimage = rectified_lf.subaperture(u,v)
                aperture_avg[v, u] = subimage.mean()

        if aperture_avg.max() == 0:
            print 'ERROR: radiometry calibration failed.  No non-zero pixels found in the calibrated radiometry image!'
            exit(1)

        # OLD RADIOMETRY (still used in ray optics model... though deprecated and we should remove!)
        # aperture_avg /= aperture_avg.max() # Make sure radiometry is scaled <= 1.0
        self.radiometry = aperture_avg

        # "FULL" RADIOMETRY: Create a full radiometric correction for every pixel in the LF.
        #
        # If the user has supplied a radiometric calibration image, we
        # use it to correct radiometry over both angles and space in
        # the light field.
        if radiometry_frame_file != None:
            print '\t    calibrating using full radiometric model.'
            self.rectified_radiometry = rectified_lf.asimage(representation = LightField.TILED_LENSLET)

        else:
            print '\t    calibrating using basic radiometric model.'
            self.rectified_radiometry = np.tile(aperture_avg, (self.nt, self.ns))

        self.radiometric_correction = self.rectified_radiometry / (self.ideal_lf + 1e-16)
        self.radiometric_correction[np.nonzero(self.ideal_lf == 0)] = 0.0    # Prevent NaNs!
        self.radiometric_correction[np.nonzero(self.rectified_radiometry == 0)] = 0.0   # Prevent NaNs!

        # Sometimes there are edge artifacts due to a strip of black
        # pixels on the left and right sides of the Neo frames.  This
        # prevents us from artificially enhancing noise in these
        # regions.
        self.radiometric_correction[ 0:self.nu,:] = 0.0;
        self.radiometric_correction[ -self.nu:,:] = 0.0;
        self.radiometric_correction[ :, 0:self.nu] = 0.0;
        self.radiometric_correction[ :, -self.nu:] = 0.0;
        
        #save_image("ideal.tif", self.ideal_lf.astype(np.float32))
        #save_image("actual.tif", self.rectified_radiometryxo.astype(np.float32))
        #save_image("correction.tif", self.radiometric_correction)

        # Set the flag that radiometry is now calibrated
        self.radiometry_is_calibrated = True
    
    # ----------------------- RAYSPREAD DATABASE --------------------------------
    
    def generate_raydb(self, num_slices, um_per_slice, supersample_factor = 1, z_center = 0.0,
                       num_threads = 1, use_geometric_optics = True, voxels_as_points = False):
        from lflib.optics import LensletArray

        self.lenslet_array = LensletArray(self.nu, self.nv, self.ns, self.nt,
                                          self.ulens_pitch, self.pixel_size,
                                          self.ulens_focal_length, self.ulens_focal_distance,
                                          self.objective_magnification, self.objective_na, self.medium_index,
                                          self.tubelens_focal_length,
                                          self.ulens_fill_factor, self.pixel_fill_factor,
                                          self.ulens_profile, self.center_wavelength)
        if use_geometric_optics:
            self.rayspread_db = self.lenslet_array.rayspread_database(num_slices, z_center,
                                                                 um_per_slice, self.radiometry,
                                                                 supersample_factor, num_threads)
        else:
             self.psf_db = self.lenslet_array.physical_rayspread_database(num_slices, z_center,
                                                                          um_per_slice, self.radiometry,
                                                                          supersample_factor,
                                                                          voxels_as_points, num_threads)
        
        # filling up the volume calibration data
        self.num_z_slices   = num_slices
        self.um_per_z_slice = um_per_slice
        self.z_center       = z_center
        
        # Set the flag that the rayspread db is now generated
        self.rayspread_db_generated = True

