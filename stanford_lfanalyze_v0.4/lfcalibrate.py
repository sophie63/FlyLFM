#!/usr/bin/env python2.7

# __BEGIN_LICENSE__
#
# Copyright 2012 Stanford University.  All rights reserved.
#
# __END_LICENSE__

# calibrate.py
#
# Usage: calibrate.py <calibration_dir>
#
# This script re-runs the calibration from a set of calibration images
# captured by the uScope GUI.  This script is mostly useful for
# debugging the calibration routines offline (i.e. away from the
# microscope), but can be used to update old calibration files from
# early experiments as well.
import sys, os, math
from lflib.imageio import load_image, save_image
import numpy as np
import glob

def avg_images(image_files):
    """
    Averages a set of images passed.
    Numerical error will become an issue with large number of files.
    """
    if len(image_files) == 0:
        return None
    im = load_image(image_files[0])
    im_type = im.dtype
    im = im.astype('float32')
    # numerical error will become a problem with large numbers of files
    for ifile in image_files[1:]:
        im = im + load_image(ifile, dtype='float32')
    return np.round(im/len(image_files)).astype(im_type)

if __name__ == "__main__":
    import lflib
    print 'LFcalibrate v%s' % (lflib.version)

    # Parse command line options
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-o", "--output-filename", dest="output_filename",
                      help="Specify the name of the calibration file.")
    parser.add_option('', "--synthetic",
                      action="store_true", dest="synthetic_lf", default=False,
                      help="Use this option to create a synthetic light field (i.e. with no calibration image")

    parser.add_option('', "--use-ray-optics",
                      action="store_true", dest="use_ray_optics", default=False,
                      help="Use the less accurate ray optics model rather than wave optics model.")
    parser.add_option('', "--voxels-as-points",
                      action="store_true", dest="voxels_as_points", default=False,
                      help="Treat each voxel as an ideal point source.  This turns of numerical integration that gives the voxel spatial extent (which can be important for anti-aliasing.")

    # Calibration routine parameters
    parser.add_option("", "--dark-frame", dest="dark_frame_file", default=None,
                      help="Specify the a dark frame image to subtract from the input light-field before processing.  (This makes radiometric calibration more accurate.)")

    parser.add_option("", "--radiometry-frame", dest="radiometry_frame_file", default=None,
                      help="Specify the a radiometry frame to use for radiometric correction.  If no frame is specified, then no radiometric correction is carried out.")
    parser.add_option("", "--align-radiometry", action="store_true", dest="align_radiometry", default=False,
                      help="Align the radiometry image automatically to the geometric calibration image.  (use this option when the radiometry frame has been \"bumped\" before imaging begins.")

    # Optical parameters
    parser.add_option('', "--pitch", dest="ulens_pitch", type=float, default=None,
                      help="Specify the microlens pitch (in microns).")
    parser.add_option('', "--pixel-size", dest="pixel_size", type=float, default=None,
                      help="Specify the size of a pixel on the sensor, taking magnification due to relay optics into account (in microns).")
    parser.add_option('', "--focal-length", dest="ulens_focal_length", type=float, default=2432.72,
                      help="Specify the microlens focal length (in microns).")
    parser.add_option('', "--ulens-focal-distance", dest="ulens_focal_distance", type=float, default=None,
                      help="Specify the microlens focal distance (in microns).  If you do not specify a value, it is assumed that the focal distance is equal to the focal length.")
    parser.add_option('', "--magnification", dest="objective_magnification", type=int, default=20,
                      help="Specify the objective magnification.")
    parser.add_option('', "--na", dest="objective_na", type=float, default = 0.5,
                      help="Specify the objective numerical aperture.")
    parser.add_option('', "--tubelens-focal-length", dest="tubelens_focal_length", type=float, default = 200.0,
                      help="Tube lens focal length (in millimeters).")
    parser.add_option('', "--wavelength", dest="center_wavelength", type=float, default = 509,
                      help="Center wavelength of emission spectrum of the sample (nm).")
    parser.add_option('', "--medium-index", dest="medium_index", type=float, default = 1.33,
                      help="Set the index of refraction of the medium.")
    parser.add_option('', "--ulens-fill-factor", dest="ulens_fill_factor", type=float, default=1.0,
                      help="Specify the microlens fill factor (e.g. 1.0, 0.7, ...).")
    parser.add_option('', "--pixel-fill-factor", dest="pixel_fill_factor", type=float, default=1.0,
                      help="Specify the pixel fill factor (e.g. 1.0, 0.7, ...).")
    parser.add_option('', "--ulens-profile", dest="ulens_profile", default='rect',
                      help="Specify the shape of the microlens apertures.  Options include: ['rect', 'circ']")


    # Volume parameters
    parser.add_option('', "--num-slices", dest="num_slices", type=int, default=30,
                      help="Set the number of slices to produce in the output stacks.")
    parser.add_option('', "--um-per-slice", dest="um_per_slice", type=float, default=10.0,
                      help="Set the thickness of each slice (in um).")
    parser.add_option('', "--z-center", dest="z_center", type=float, default=0.0,
                      help="Set the offset for the central z slice (in um).")
    parser.add_option('', "--supersample", dest="supersample", type=int, default= 1,
                      help="Supersample the light field volume.  This results in a higher resolution reconstruction up to a point, and interpolation after that point.")

    # Geometric calibration Options
    parser.add_option("", "--affine-alignment", action="store_true", dest="affine_alignment", default=False,
                      help="Use affine warp for correcting geometric distortion (default is cubic).")
    parser.add_option("", "--isometry-alignment", action="store_true", dest="isometry_alignment", default=False,
                      help="Use isometry warp for correcting geometric distortion (default is cubic).")
    parser.add_option("--chief-ray", action="store_true", dest="chief_ray_image", default=False,
                      help="Use this flag to indicate that the calibration frame is a chief ray image.")

    # Synthetic parameters
    parser.add_option('', "--ns", dest="ns", type=int, default=50,
                      help="Set the lenslets in s direction.")
    parser.add_option('', "--nt", dest="nt", type=int, default=50,
                      help="Set the lenslets in t direction.")
    
    # Other Options
    parser.add_option('', "--crop-center-lenslets",
                      action="store_true", dest="crop_center_lenslets", default=False,
                      help="For severe aperture vignetting (high NA objectives), use only center lenslets for calibration, and extrapolate outwards.")
    parser.add_option('', "--skip-alignment",
                      action="store_true", dest="skip_alignment", default=False,
                      help="Skip the alignment step during geometric calibration (useful if you are working with an already-rectified light field or a synthetic light field.")
    parser.add_option("", "--skip-subpixel-alignment", action="store_true", dest="skip_subpixel_alignment", default=False,
                      help="Skip subpixel alignment for determining lenslet centers.")
    parser.add_option('', "--num-threads", dest="num_threads", type=int, default=10,
                      help="Set the number of CPU threads to use when generating the raydb.")
    parser.add_option("", "--pinhole", dest="pinhole_filename", default=None,
                      help="After calibrating, save the rectified light field as a rectified sub-aperture image.")
    parser.add_option("", "--lenslet", dest="lenslet_filename", default=None,
                      help="After calibrating, save the rectified light field as a rectified lenslet image.")
    parser.add_option('-d', "--debug",
                      action="store_true", dest="debug", default=False,
                      help="Save debug images.")

    (options, args) = parser.parse_args()

    # If no focal distance is supplied, then set it (by default) to be equal to the ulens focal length.
    if options.ulens_focal_distance == None:
        options.ulens_focal_distance = options.ulens_focal_length

    if not options.synthetic_lf and len(args) != 1:
        print 'You must supply exactly one calibration image.\n'
        sys.exit(1)
    calibration_filename = args[0]

    if options.pixel_size == None or options.ulens_pitch == None:
        print 'Please supply necessary pixel per lenslet information via the \'--pitch\' and \'--pixel-size\' options.'
        sys.exit(1)

    if options.synthetic_lf:
        ns = options.ns
        nt = options.nt
        nu = nv = int(np.ceil(float(options.ulens_pitch) / options.pixel_size))
        synthetic_lf = 65535 * np.ones((nt*nv, ns*nu), dtype=np.uint16)
        save_image(calibration_filename, synthetic_lf, dtype=np.uint16)
        options.skip_alignment = True

    # Default output filename has a -RECTIFIED suffix
    if not options.output_filename:
        fileName, fileExtension = os.path.splitext(calibration_filename)
        output_filename = fileName + '.lfc'
    else:
        output_filename = options.output_filename

    # Check if dark-frame or radiometry-frame are regular expressions referring to multiple files.
    # If so, save an average image as the dark/radiometry frame
    if options.dark_frame_file is not None and len(glob.glob(options.dark_frame_file)) > 1:
        dark_frame_files = glob.glob(options.dark_frame_file)
        avg_dark_frame = avg_images(dark_frame_files)
        options.dark_frame_file = os.path.dirname(output_filename) + os.sep + 'darkframe_avg' + os.path.splitext(dark_frame_files[0])[1]
        save_image(options.dark_frame_file, avg_dark_frame)

    if options.radiometry_frame_file is not None and len(glob.glob(options.radiometry_frame_file)) > 1:
        radiometry_frame_files = glob.glob(options.radiometry_frame_file)
        avg_radiometry_frame = avg_images(radiometry_frame_files)
        options.radiometry_frame_file = os.path.dirname(output_filename) + os.sep + 'radiometryframe_avg' + os.path.splitext(radiometry_frame_files[0])[1]
        save_image(options.radiometry_frame_file, avg_radiometry_frame)

    # Create a new calibration object
    from lflib.calibration import LightFieldCalibration
    

    # FOR DEBUGGING: Load a previous calibration
    #
    #lfcal = LightFieldCalibration.load(output_filename)

    from lflib.calibration.imaging import CalibrationAlignmentMethods
    if options.affine_alignment:
        calibration_method = CalibrationAlignmentMethods.CALIBRATION_AFFINE_ALIGNMENT
    elif options.isometry_alignment:
        calibration_method = CalibrationAlignmentMethods.CALIBRATION_ISOMETRY_ALIGNMENT
    else:
        calibration_method = CalibrationAlignmentMethods.CALIBRATION_CUBIC_ALIGNMENT

    lfcal = LightFieldCalibration(options.ulens_focal_length, options.ulens_focal_distance,
                                  options.ulens_pitch, options.pixel_size,
                                  options.objective_magnification, options.objective_na, options.medium_index,
                                  options.tubelens_focal_length,
                                  options.ulens_fill_factor, options.pixel_fill_factor,
                                  options.ulens_profile, options.center_wavelength,
                                  calibration_method)

    # STEP 1 : MEASURE THE GEOMETRIC DISTORTION
    # 
    # This routine computes an affine transform that squares the
    # lenslet array to a nice, regularly space grid.
    lfcal.calibrate_geometry(calibration_filename,
                             skip_alignment = options.skip_alignment,
                             skip_subpixel_alignment = options.skip_subpixel_alignment,
                             debug_mode = options.debug,
                             chief_ray_image = options.chief_ray_image,
                             radiometry_file = options.radiometry_frame_file,

                             align_radiometry = options.align_radiometry,
                             crop_center_lenslets = options.crop_center_lenslets)
    print ('   Calibrated light field has [ %d x %d ] ray samples and [ %d x %d ] spatial samples.' %
           (lfcal.nu, lfcal.nv, lfcal.ns, lfcal.nt))

    # Optionally, create a rectified sub-aperture image
    if (options.pinhole_filename):
        from lflib.lightfield import LightField
        im = load_image(calibration_filename, dtype=np.float32, normalize = False)
        lf = lfcal.rectify_lf(im).asimage(LightField.TILED_SUBAPERTURE)
        save_image(options.pinhole_filename, lf/lf.max() * 65535, dtype=np.uint16)

    # Optionally, create a rectified lenslet image
    if (options.lenslet_filename):
        from lflib.lightfield import LightField
        im = load_image(calibration_filename, dtype=np.float32, normalize = False)
        lf = lfcal.rectify_lf(im).asimage(LightField.TILED_LENSLET)
        save_image(options.lenslet_filename, lf/lf.max() * 65535, dtype=np.uint16)

    # For debugging
    #print "DEBUG MODE ON!!!!!"
    #raise SystemExit


    # STEP 2 : Compute rayspread database
    #
    # The rayspread database is a look-up table that serves as the
    # optical model of forward and back-projection of the light field.
    print '-> Generating light field psf database.  (This may take a little while...)'
    lfcal.generate_raydb(options.num_slices, options.um_per_slice,
                         options.supersample, options.z_center, options.num_threads,
                         use_geometric_optics = options.use_ray_optics,
                         voxels_as_points = options.voxels_as_points)

    # STEP 3 : MEASURE THE APERTURE PLANE VIGNETTING FOR THIS LENSLET IMAGE
    #
    # The vignetting function can be used for deconvolution.
    
    # First we determine a reasonable number of pixels per lenslet to
    # use.  This must be the same scheme used in lfstack.py and
    # elsewhere.  It's a little dangerous here to be accessing
    # coefficient directly... we should veil this in some layer of
    # abstraction soon!
    print '-> Calibrating radiometry using ', options.radiometry_frame_file
    from lflib.lfexceptions import ZeroImageException
    try:
        lfcal.calibrate_radiometry(calibration_filename, 
                radiometry_frame_file = options.radiometry_frame_file,
                dark_frame_file = options.dark_frame_file)
        # Save the result
        lfcal.save(output_filename);
        lfcal.print_summary()
        
    except ZeroImageException:
        print "ERROR: calibrating against a blank radiometric image"
