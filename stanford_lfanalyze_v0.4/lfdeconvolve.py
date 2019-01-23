#!/usr/bin/env python2.7

# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2013 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

# lflib imports 
import lflib
from lflib.imageio import load_image, save_image
from lflib.lightfield import LightField
from lflib.calibration import LightFieldCalibration
from lflib.util import ensure_path

# major libraries
import numpy as np
import h5py
import math, sys
import cv
import os
import subprocess
import time
import tempfile
import socket
import argparse

#----------------------------------------------------------------------------------

class SCPFailure(Exception):
    def __init__(self, filepath, error_output):
        Exception.__init__(self, filepath, error_output)
        self.filepath = filepath
        self.errout = error_output
    def __str__(self):
        return """
        #####################
        SCP failed for %s.
        Process output was:
        %s
        #####################
        """ % (self.filepath, self.errout)

def copy_file(source_path, dest_path):
    "Copy a file via SCP, retrying a few times if necessary."
    retries = 5
    retry_interval_secs = 1
    while retries > 0:
        tic = time.time()
        copyproc = subprocess.Popen(['scp', '-o StrictHostKeyChecking=no', '-i', args.private_fn, source_path, dest_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = copyproc.communicate()
        if copyproc.returncode == 0:
            print '\t\t--> file tranfer completed in %f secs' % (time.time() - tic)
            break
        time.sleep(retry_interval_secs)
        retry_interval_secs *= 2 # backoff retry time
        print '\t\t--> file transfer failed.  Retrying.'
    else:
        # retries exhausted.  give up.
        raise SCPFailure(calibration_filename, err)


def retrieve_calibration_file(calibration_filename, id=''):

    if ':' in calibration_filename and '@' in calibration_filename:
        print '\t--> fetching calibration file on remote host...'
        cache_filename = tempfile.gettempdir() + os.sep + 'tmp_lfc_cache_' + str(hash(calibration_filename)) + '_' + id + '.lfc'

        # If the file is already present on the filesystem, then we can use it immediately.
        if os.path.isfile(cache_filename):
            print '\t\t--> using locally cached copy %s.'%cache_filename
            calibration_filename = cache_filename

        # If the file is not present, then it needs to be downloaded.
        # Only one worker can do this at a time.  The rest must wait
        # for the file to be downloaded.  We use a dummy file on the
        # filesystem as a concurrency "lock" to ensure that only one
        # worker downloads the file.
        else:
            try:
                # Try to open a file to use as a concurrency lock.  Only the first process to try will succeed.  The rest will fail.
                fd_concurrency_lock = os.open(cache_filename + ".lock", os.O_CREAT|os.O_EXCL|os.O_RDWR)

                # If the worker has gotten this far, it has the concurrency lock.
                print '\t\t--> using scp. temporary transfer location: %s.'%(cache_filename+'.tmp')
                copy_file(calibration_filename, cache_filename+'.tmp')

                # Rename the file
                os.rename(cache_filename+'.tmp', cache_filename)

                # Remove the concurrency lock.  This should signal to
                # other processes that they can safely use the
                # downloaded lfc file.
                os.close(fd_concurrency_lock)
                os.remove(cache_filename + ".lock")
                print '\t\t--> transfer complete. Cached at: %s.' % cache_filename
                
            except OSError:
                print '\t\t--> another process appears to be retriving the file.  Waiting for transfer to complete.'
                numAttempts = 0
                while os.path.isfile(cache_filename + ".lock") and numAttempts < 10000:
                    time.sleep(5)
                    numAttempts+=1
                    if numAttempts>1000:
                        raise Exception('completion of transfer of calibration file by another process seems to stalled.')
            
            calibration_filename = cache_filename

    return calibration_filename

def do_deconvolve(args):

    filename = args.input_file

    # Default calibration filename has a *.lfc suffix, and the same prefix
    if not args.calibration_file:
        fileName, fileExtension = os.path.splitext(filename)
        calibration_file = fileName + '.lfc'
    else:
        calibration_file = args.calibration_file

    # Default output filename has a -RECTIFIED suffix
    if not args.output_filename:
        fileName, fileExtension = os.path.splitext(filename)
        output_filename = fileName + '-STACK.tif'
    else:
        output_filename = args.output_filename

    print '\t--> hostname:{host}'.format(host=socket.gethostname())
    print '\t--> specified gpu-id:{gpuid}'.format(gpuid=args.gpu_id)

    #check if output filename appears to be formatted as user@host:path.                                                                                                            
    remote_output_fn = None
    if ':' in output_filename and '@' in output_filename:
        remote_output_fn = output_filename
        output_filename = tempfile.gettempdir() + os.sep + 'deconv_tmpout_' + str(hash(filename)) + os.path.splitext(output_filename)[1]
        print '\t--> output file location is on remote host. Output will temporarily appear at %s.'%output_filename

    # Loadim the calibration data 
    calibration_file = retrieve_calibration_file(calibration_file, id=str(args.gpu_id))
    lfcal = LightFieldCalibration.load(calibration_file) 
    print '\t--> loaded calibration file.'

    # Load the raw input data
    if ':' in filename and '@' in filename:
        tmpname = tempfile.gettempdir() + os.sep + 'deconv_tmpin_' + str(hash(filename)) + os.path.splitext(filename)[1]
        print '\t--> input file is on remote host. Transferring temporarily to  %s.'%tmpname 
        copy_file(filename, tmpname)
        im = load_image(tmpname, dtype=np.float32, normalize=False)
        os.remove(tmpname)
    else:
        im = load_image(filename, dtype=np.float32, normalize = False)
    print '\t--> %s opened.  Pixel values range: [%d, %d]' % (filename, int(im.min()), int(im.max()))

    # Perform dark frame subtraction
    from lflib.lfexceptions import ZeroImageException
    try:
        im = lfcal.subtract_dark_frame(im)
        print '\t    Dark frame subtracted.  Pixel values range: [%f, %f]' % (im.min(), im.max())
        lf_zeros = False
    except ZeroImageException:
        print "\t    A frame with no light pixels was found, but it's no big deal"
        lf_zeros = True
    
    # Rectify the image
    lf = lfcal.rectify_lf(im)

#    save_image('/home/broxton/debug_retified_lenslet.tif', lf.asimage(LightField.TILED_LENSLET), dtype=np.float32)
#    save_image('/home/broxton/debug_retified_subap.tif', lf.asimage(LightField.TILED_SUBAPERTURE), dtype=np.float32)
#    save_image('/home/broxton/debug_retified_radiometry.tif', lfcal.rectified_radiometry, dtype=np.float32)

    # Save a little bit of verboseness in the code below by extracting the appropriate 'db' object.
    if lfcal.psf_db != None:
        print '\t    Using wave optic psf db.'
        db = lfcal.psf_db
    else:
        print '\t    Using rayspread db.'
        db = lfcal.rayspread_db
    
    # Initialize the volume as a plain focal stack.  We normalize by the weights as well.
    from lflib.volume import LightFieldProjection
    lfproj = LightFieldProjection(lfcal.rayspread_db, lfcal.psf_db,
                                  disable_gpu = args.disable_gpu, gpu_id = args.gpu_id)

    # Enable radiometry correction
    lfproj.set_premultiplier(lfcal.radiometric_correction)

    # DEBUG - MANUAL RADIOMETRY FOR TRYING VARIOUS STUFF OUT
    # lf_ones = lfproj.project(db.ones_volume())
    # ideal_lf = lf_ones.asimage(representation = LightField.TILED_LENSLET)
    # print ideal_lf.min(), ideal_lf.max()
    # rectified_radiometry = lfcal.rectified_radiometry
    # print rectified_radiometry.min(), rectified_radiometry.max()
    # radiometric_correction = rectified_radiometry / (ideal_lf + 1e-16)
    # print radiometric_correction.min(), radiometric_correction.max()
    # save_image("raddebug_ideal.tif", ideal_lf, dtype=np.float32)
    # save_image("raddebug_actual.tif", rectified_radiometry, dtype=np.float32)
    # save_image("raddebug_correction.tif", radiometric_correction, dtype=np.float32)

    # lightfield_im = lf.asimage(representation = LightField.TILED_LENSLET)


    if args.test is not None:
        print '===================================='
        if args.test == 1:
            lfproj.test_forward_backward()
        if args.test == 2:
            lfproj.test_forward_project()
        if args.test == 3:
            lfproj.test_back_project()
        print 'finished testing SIRT, quitting...'
        exit()
    
    print '-------------------------------------------------------------------'
    print 'Computing light field tomographic reconstruction of:', filename
    print 'conv_thresh:', args.conv_thresh
    print '   max_iter:', args.max_iter
    print '     lambda: ', args.regularization_lambda

    # All of the solvers below take as arguments a system matrix
    # operator 'A' (implemented as a numpy LinearOperator) and a 'b'
    # matrix, which is the light field sub-aperture image transformed
    # into a vector.  Some routines also take a starting volume vector
    # 'x'.
    #
    # We take the resulting vector 'x' and transform it back into a
    # volume below.
    #
    from lflib.linear_operators import LightFieldOperator
    nrays = db.ns*db.nu*db.nt*db.nv
    nvoxels = db.nx*db.ny*db.nz
    A_lfop = LightFieldOperator(lfproj, db)
    A_operator = A_lfop.as_linear_operator(nrays, nvoxels)

    # Trim out entries in the light field that we are ignoring because
    # they are too close to the edge of the NA of the lenslet.  This
    # make our residuals more accurate below.
    lf = lfcal.mask_trimmed_angles(lf)

    # Generate the b vector, which contains the observed lightfield;
    # and the initial volume x containing all zeros.
    im_subaperture = lf.asimage(representation = LightField.TILED_SUBAPERTURE)
    b_vec = np.reshape(im_subaperture, np.prod(im_subaperture.shape))

    # DEBUG  -  TEST OUT INTERCEPT IDEA
    #
    # This code seems to work well for removing some of the edge
    # artifacts in our volumes, but I would still classify it as
    # "beta" code for now.  I'm leaving it here but we should consider
    # whether this is exactly what we want to be doing.
    #
    TEST_INTERCEPT_CODE = False
    if TEST_INTERCEPT_CODE:
        lf_ones = lfproj.project(db.ones_volume()).asimage(LightField.TILED_SUBAPERTURE)
        lf_ones /= lf_ones.max()

        from lflib.linear_operators import LightFieldOperatorWithFullIntercept
        A_lfop = LightFieldOperatorWithFullIntercept(lfproj, db)
        A_operator = A_lfop.as_linear_operator()

        # from lflib.linear_operators import LightFieldOperatorWithIntercept
        # A_lfop = LightFieldOperatorWithIntercept(lfproj, db, np.reshape(lf_ones, np.prod(lf_ones.shape)))
        # A_operator = A_lfop.as_linear_operator()
    # /DEBUG

    if (args.benchmark):
        lfproj.compare_cpu_gpu_performance()
        raise SystemExit
    
    if lf_zeros:
        vol = lfcal.psf_db.empty_volume()
        x_vec = np.reshape(vol, np.prod(vol.shape)).astype(np.float32)

    # If the user has requested a focal stack, we stop here and return
    # the current volume.
    elif args.focalstack:
        print 'Computing focal stack reconstruction of:', filename

        # Weight come from one forward projection followed by one back projection.
        vol_weights = lfproj.backproject(lfproj.project(db.ones_volume()))

        # This multiplier removes the grid artifact entirely from the
        # volume.  This only seems to work for the focal stack,
        # though.
        pm = 1.0 / vol_weights
        pm[np.nonzero(vol_weights == 0)] = 0.0               # Prevent NaNs!
        lfproj.set_postmultiplier(pm);

        vol = lfproj.backproject(lf)
        TEST_INTERCEPT_CODE = False                          # Need to disable this for focal stack.
        x_vec = np.reshape(vol, np.prod(vol.shape)).astype(np.float32)
        
    elif args.solver == 'amp':
        from lflib.solvers.amp import amp_reconstruction
        if args.conv_thresh == 0.0: args.conv_thresh = 1e-10
        args.alpha = 1.0
        vol, multiscale_coefs = amp_reconstruction(lfcal, lf, args.alpha,
                                                   args.conv_thresh, args.max_iter,
                                                   args.regularization_lambda,
                                                   delta = 1.0, 
                                                   debug = args.debug,
                                                   disable_gpu = args.disable_gpu,
                                                   gpu_id = args.gpu_id,
                                   debug_path = args.output_filename.split('.tif')[0] + '_',
                                    multiscale_smoothing=args.multiscale_smoothing)
        x_vec = np.reshape(vol, np.prod(vol.shape))

    elif args.solver == 'admm_huber':
        print "WARNING: ADMM Huber solver is experimental!"
        from lflib.solvers.admm import admm_huber_reconstruction
        if args.conv_thresh == 0: convergence_threshold = 1e-4
        vol = admm_huber_reconstruction(lfcal, lf, args.alpha,
                                        args.conv_thresh, args.max_iter,
                                        args.regularization_lambda,
                                        disable_gpu = args.disable_gpu,
                                        gpu_id = args.gpu_id)
        x_vec = np.reshape(vol, np.prod(vol.shape))

    elif args.solver == 'admm_tv':
        print "WARNING: ADMM Total Variation solver is experimental!"
        from lflib.solvers.admm import admm_total_variation_reconstruction
        if args.conv_thresh == 0: convergence_threshold = 1e-4
        vol = admm_total_variation_reconstruction(lfcal, lf, args.alpha,
                                                  args.conv_thresh, args.max_iter,
                                                  args.regularization_lambda,
                                                  args.regularization_lambda2,
                                                  disable_gpu = args.disable_gpu,
                                                  gpu_id = args.gpu_id)
        x_vec = np.reshape(vol, np.prod(vol.shape))

    elif args.solver == 'cg':
        from lflib.solvers.conjugate_gradient import conjugate_gradient_reconstruction
        if args.conv_thresh == 0: convergence_threshold = 1e-9
        vol = conjugate_gradient_reconstruction(lfcal, lf, 
                                                args.conv_thresh, args.max_iter,
                                                args.regularization_lambda,
                                                disable_gpu = args.disable_gpu,
                                                gpu_id = args.gpu_id)
        x_vec = np.reshape(vol, np.prod(vol.shape))

    elif args.solver == 'rl':
        from lflib.solvers.richardson_lucy import richardson_lucy_reconstruction
        if args.conv_thresh == 0: args.conv_thresh = 1e-6
        x_vec = richardson_lucy_reconstruction(A_operator, b_vec,
                                               Rtol = args.conv_thresh,
                                               max_iter = args.max_iter,
                                               beta = args.background_level,
                                               sigmaSq = args.readnoise_variance)

    elif args.solver == 'rl_multiscale':
        from lflib.solvers.richardson_lucy_multiscale import richardson_lucy_multiscale_reconstruction
        vol_shape = (db.ny, db.nx, db.nz)
        if args.conv_thresh == 0: args.conv_thresh = 1e-6
        x_vec = richardson_lucy_multiscale_reconstruction(A_operator, b_vec, vol_shape,
                                                          Rtol = args.conv_thresh,
                                                          max_iter = args.max_iter,
                                                          beta = args.background_level,
                                                          sigmaSq = args.readnoise_variance)


    elif args.solver == 'mrnsd':
        from lflib.solvers.mrnsd import wmrnsd_reconstruction
        if args.conv_thresh == 0: convergence_threshold = 1e-6
        x_vec = wmrnsd_reconstruction(A_operator, b_vec,
                                     Rtol = convergence_threshold,
                                     max_iter = args.max_iter,
                                      beta = args.background_level,
                                      sigmaSq = args.readnoise_variance)


    elif args.solver == 'direct':
        print "WARNING: Direct deconvolution will only work on smaller problems!"
        from lflib.solvers.direct import direct_reconstruction
        vol = direct_reconstruction(ray_db, args. lfcal.lenslet_array, aperture_radiometry, lfproj)
        x_vec = np.reshape(vol, np.prod(vol.shape))

    elif args.solver == 'lsqr':
        from lflib.solvers.lsqr import lsqr_reconstruction
        if args.conv_thresh == 0: convergence_threshold = 1e-9
        x_vec = lsqr_reconstruction(A_operator, b_vec,
                                    args.conv_thresh, args.max_iter,
                                    args.regularization_lambda)

    elif args.solver == 'kspace':
        print "WARNING: Use for testing only!"
        # use for testing only!
        # um_per_slice is not actually used right now; 
        # a volume with arbitrary z-spacing is returned
        from  lflib.kspace_deconvolution import lfdeconvolve_kspace
        apertureimage_lf = lf.asimage(representation = LightField.TILED_SUBAPERTURE)
        vol = lfdeconvolve_kspace( nt, ns, ray_db.nz, 0, nu, nv,
                                   lfcal.magnification,
                                   lfcal.na, apertureimage_lf, mediumRI = lfcal.medium_index,
                                   umPerLenslet = lfcal.array_pitch, 
                                   focalstack = args.focalstack  )
        x_vec = np.reshape(vol, np.prod(vol.shape))

    elif args.solver == 'sirt':
        from lflib.solvers.sirt import sirt_reconstruction
        if args.conv_thresh == 0: convergence_threshold = 5e-5
        (x_vec, residuals) = sirt_reconstruction(A_operator, b_vec,
                                                 args.alpha, args.conv_thresh, args.max_iter)
    
    else:
        raise ValueError("The reconstruction method specified is not an option, " + 
                         "available methods are: "
                         " Approximate Message Passing (with optional multiscale denoising ('amp'), " +
                         " Alternating Direction Method of Multipliers with Huber loss ('admm_huber'), "+
                         " Alternating Direction Method of Multipliers with TV penalty ('admm_tv), " +
                         " Conjugate Gradient ('cg'), " + 
                         " Direct method with Cholesky factorization ('direct'), " +
                         " Least Squares QR ('lsqr'), " + 
                         " K-space deconvolution ('kspace'), " +
                         " Simultaneous Iterative Reconstruction Technique ('sirt'), " +
                         " MRNSD ('mrnsd')," + 
                         " and Richardson-Lucy ('rl').")


    # ============================ RESHAPE & CLEAN UP VOLUME ==============================

    # DEBUG  -  TEST OUT INTERCEPT IDEA
    if TEST_INTERCEPT_CODE:
        background_term = x_vec[db.nvoxels:]
        background_im = np.reshape(background_term, (db.nt*db.nv, db.ns*db.nu))
        # save_image("debug_background.tif", background_im, dtype=np.float32)

        x_vec = x_vec[0:db.nvoxels]
    # /DEBUG

    vol = np.reshape(x_vec, (db.ny, db.nx, db.nz))
    
    # Slight hack: zero out the outermost XY "shell" of pixels, since
    # these are often subject to radiometry artifacts.
    min_val = vol[db.supersample_factor:-db.supersample_factor,
                  db.supersample_factor:-db.supersample_factor, :].min()
    print '\t--> Replacing border values with min value: ', min_val
    vol[0:db.supersample_factor, :, :] = min_val
    vol[-db.supersample_factor:, :, :] = min_val
    vol[:, 0:db.supersample_factor, :] = min_val
    vol[:, -db.supersample_factor:, :] = min_val

    # ================================= SAVE RESULTS ==================================

    # (optionally) remove grid artifacts
    if args.remove_grid:
        from lflib.postfilter import remove_light_field_grid
        if lfcal.psf_db != None:
            print '\t--> Removing grid artifacts in light field image using spectral median filter.'
            vol = remove_light_field_grid(vol, lfcal.psf_db.supersample_factor)
        else:
            print '\t--> Skipping grid artifact removal.  Not needed for ray optics reconstructions.'

    # Create a new hdf5 file for the output
    if os.path.splitext(output_filename)[1].lower() == ".h5":
        f = h5py.File(output_filename,'w')
        f.create_dataset('timeseries_volume', data=np.reshape(vol, np.prod(vol.shape), order='f'))
        vol_shape = np.array([vol.shape[2], vol.shape[1], vol.shape[0]])
        f.create_dataset('vol_shape', data=vol_shape)
        print '    Saving', f
        f.close()

    # Or, if the file suffix is *.tif, then we save a tiff stack
    else:
        # Save volume
        save_image(output_filename, vol)

        # If multiscale decomposition is used, save multiscale images
        if args.multiscale_smoothing and args.save_multiscale:
            print "\t--> Saving multiscale time series volumes as TIFF stacks."
            for i in xrange(len(multiscale_coefs)):
                base_path = os.path.split(output_filename)[0]
                base_name = os.path.splitext( os.path.split(output_filename)[1])[0]
                out_path = os.path.join(base_path, "multiscale_stack_cache_"+str(i))

                # Create output directory if if does not already exist
                if not os.path.isdir(out_path):
                    print('\t--> Creating multiscale stack directory: '+out_path)
                    os.makedirs(out_path)

                multiscale_output_filename = os.path.join(out_path, base_name+'-MULTISCALE'+
                                                          str(i)+'.tif') 
                save_image(multiscale_output_filename, multiscale_coefs[i])

    # If necessary, transfer output file to remote host and delete local copy.                                                                                                      
    if remote_output_fn:      
        copy_file(output_filename, remote_output_fn)
        os.remove(output_filename)    

    # Optionally, create a "filtered" pinhole image 
    if (args.pinhole_filename):
        lf = lfproj.project(vol).asimage(LightField.TILED_SUBAPERTURE)
        save_image(args.pinhole_filename, lf/lf.max() * 65535, dtype=np.uint16)

#----------------------------------------------------------------------------------
#                                        MAIN
#----------------------------------------------------------------------------------
if __name__ == "__main__":
    print 'LFdeconvolve v%s' % (lflib.version)

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', 
                        help="You must supply at least one light field image "
                        "to deconvolve.")
    parser.add_argument("-o", "--output-file", dest="output_filename",
                        help="Specify the output filename.")
    parser.add_argument("-c", "--calibration-file", dest="calibration_file", 
                        help="Specify the calibration file to use for rectification.")
    key_def = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                           'enlightenment_c3')
    parser.add_argument("--private-key", dest="private_fn", type=str, 
                        default=key_def, 
                        help="Specify the private key file for remote transfers.")
    parser.add_argument("--cov-directory", dest="cov_directory",
                        help="Specify the directory where ADMM covariance matrices are saved.")

    # Algorithm selection
    parser.add_argument("--solver", dest="solver", default='rl',
                      help= "Available reconstruction methods are: "
                         " Approximate Message Passing (with optional multiscale denoising ('amp'), " +
                         " Alternating Direction Method of Multipliers with Huber loss ('admm_huber'), "+
                         " Alternating Direction Method of Multipliers with TV penalty ('admm_tv), " +
                         " Conjugate Gradient ('cg'), " + 
                         " Direct method with Cholesky factorization ('direct'), " +
                         " Least Squares QR ('lsqr'), " + 
                         " K-space deconvolution ('kspace'), " +
                         " Simultaneous Iterative Reconstruction Technique ('sirt'), " +
                         " MRNSD ('mrnsd')," + 
                         " and Richardson-Lucy ('rl'). Default is currently 'sirt'.")

    # Algorithm-specific parameters
    parser.add_argument("--alpha", dest="alpha", type=float, default=1.6,
                        help="Relaxation parameter for SIRT-based iterative reconstruction.")
    parser.add_argument("--multiscale-smoothing", 
                        dest="multiscale_smoothing", action="store_true",
                        default=False, 
                        help="Multiscale regularization option for AMP reconstruction.")
    parser.add_argument("--save-multiscale", dest="save_multiscale", action="store_true",
                        default=False, help="Save multilevel decomposition of data.")
    # Generic parameters for iterative reconstruction routines
    parser.add_argument("--lambda", dest="regularization_lambda", 
                        type=float, default=0.0,
                        help="Regularization coefficient (behavior varies by "
                        "reconstruction algorithm)")
    parser.add_argument("--lambda2", dest="regularization_lambda2", type=float, default=0.0,
                        help=("Additional regularization coefficient. "
                            "(Behavior varies by algorithm, and not all algorithms use "
                            "two regularization coefficients.)"))
    parser.add_argument("--max-iter", dest="max_iter", type=int, default=15,
                        help="Maximum number of iterations for SIRT-based reconstruction.")
    parser.add_argument("--convergence-threshold", dest="conv_thresh", type=float, default= 0.0,
                        help=("Convergence criteria threshold, d/dt (MSE). "
                            "Try 5e-5 for SIRT, 1e-2 for TV."))

    # Noise model parameters
    parser.add_argument("--readnoise-variance", dest="readnoise_variance", type=float, default=0.0,
                      help="Set the variance of the (measured) camera read noise.")
    parser.add_argument("--background-level", dest="background_level", type=float, default=1.0,
                      help="Set the (measured) background level of the image.")

    # Assorted other parameters
    parser.add_argument("--disable-gpu",
                      action="store_true", dest="disable_gpu", default=False,
                      help="Disable GPU deconvolution, and use software implementation instead.")

    # if --gpu-id is not supplied it will default to the value of the USE_GPU_ID environment 
    # variable, if set, or 1 otherwise.
    parser.add_argument("--gpu-id", dest="gpu_id", type=int, default=int(os.environ.get('USE_GPU_ID', 1)),
                      help="Force lfdeconvolve to use a specific GPU on your system. If not supplied this will default to $USE_GPU_ID or 1")
    parser.add_argument("--focalstack",
                      action="store_true", dest="focalstack", default=False,
                      help="Turn off deconvolution and simply save a focal stack to disk.")
    parser.add_argument("--remove-grid",
                      action="store_true", dest="remove_grid", default=False,
                      help="Remove grid artifacts in light field image using spectral median filter.")
    parser.add_argument("-p", "--pinhole-file", dest="pinhole_filename", default=None,
                      help="After deconvolution, save out a deconvolved light field sub-aperture image.")
    parser.add_argument("--deconvolution-type", dest="decon_type", default='algebraic',
                      help="Choose deconvolution method. One of [algebraic, direct, admm].")
    parser.add_argument("--reg-factor", dest="reg_factor", type=float, default=100.,
                      help="Regularization parameter used in ADMM.")
    parser.add_argument("--h5py-cov-filename", dest="h5py_cov_filename", 
                      default='tests/covariance_blocks.h5',
                      help="Specify the HDF5 covariance filename.")
    parser.add_argument("--direct-type", dest="direct_type", default="covariance",
                        help="If --direct flag is set to True, specifies "
                        "whether the covariance or projection matrix method is used.")
    parser.add_argument("--benchmark",
                      action="store_true", dest="benchmark", default=False,
                      help="Compare the CPU and GPU speeds for forward & back porjection operations.")
    parser.add_argument("--test",
                      dest="test", type=int, help="Select a unit test (1-4).")
    parser.add_argument("--log-convergence",
                      action="store_true", dest="log_convergence", default=False,
                      help="For logging convergence details.")
    args = parser.parse_args()

    if args.output_filename is not None:
        ensure_path(args.output_filename)
    
    # Run the deconvolution algorithm
    do_deconvolve(args)
