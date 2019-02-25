from lflib.calibration.grid import build_grid
from lflib.calibration.ransac import ransac
from lflib.imageio import save_image

import numpy as np
import os

def difference_of_gaussians(im, radius):
    sigma=0.44248120*radius+0.01654135
    return scipy.ndimage.gaussian_filter(im, sigma) - scipy.ndimage.gaussian_filter(im, 1.6*sigma)

# Non-max Suppression Algorithm.  Returns an image with pixels that
# are not maximal within a square neighborhood zeroed out.
def nonmaxsup(im, radius, threshold):

    # Normalize the image & threshold
    im = im / im.max()
    im[np.nonzero(im<threshold)] = 0

    # Extract local maxima by performing a grey scale morphological
    # dilation and then finding points in the corner strength image that
    # match the dilated image and are also greater than the threshold.
    from scipy.ndimage import morphology
    mx = morphology.grey_dilation(im, footprint=np.ones((radius, radius)), mode='constant', cval=0)
    return im * (im >= mx);

def refine_centers_subpixel(putative_centers, lf_tophat):
    '''
    Fit a 2D parabola to the intensity values around a peak in image
    intensity.  Then return the peak of the parabola, which serves as
    relatively accurate sub-pixel estimate of the true peak location.
    '''
    # Refine lenslet centers estimates with sub-pixel accuracy.    
    pinvA = np.matrix([[ 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6],
                       [ 1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6],
                       [ 1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4],
                       [-1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,    0.0,  1.0/6,  1.0/6,  1.0/6],
                       [-1.0/6,    0.0,  1.0/6, -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0,  1.0/6],
                       [-1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9,  2.0/9, -1.0/9,  2.0/9, -1.0/9]])

    subpixel_centers = np.zeros(putative_centers.shape)
    for i in range(0, putative_centers.shape[0]):
        a0 = lf_tophat[putative_centers[i,1]-1, putative_centers[i,0]-1]
        a3 = lf_tophat[putative_centers[i,1]-1, putative_centers[i,0]+0]
        a6 = lf_tophat[putative_centers[i,1]-1, putative_centers[i,0]+1]
        a1 = lf_tophat[putative_centers[i,1]-0, putative_centers[i,0]-1]
        a4 = lf_tophat[putative_centers[i,1]-0, putative_centers[i,0]+0]
        a7 = lf_tophat[putative_centers[i,1]-0, putative_centers[i,0]+1]
        a2 = lf_tophat[putative_centers[i,1]+1, putative_centers[i,0]-1]
        a5 = lf_tophat[putative_centers[i,1]+1, putative_centers[i,0]+0]
        a8 = lf_tophat[putative_centers[i,1]+1, putative_centers[i,0]+1]
        points = np.matrix([a0, a1, a2, a3, a4, a5, a6, a7, a8]).transpose()
        c = pinvA * points
        denom = 4 * c[0] * c[1] - (c[2] * c[2])
        x = (c[2] * c[4] - 2 * c[1] * c[3]) / denom
        y = (c[2] * c[3] - 2 * c[0] * c[4]) / denom
        # print "Pixel: (%f, %f), Subpixel: (%f, %f)" % (putative_centers[i,0], putative_centers[i,1], x, y)
        subpixel_centers[i,0] = putative_centers[i,0] + x;
        subpixel_centers[i,1] = putative_centers[i,1] + y;

    return subpixel_centers


def crop_residuals_to_inside(model, residuals, calibration_image_shape):
    ul = model.eval_point([0, 0], 'f').flatten()
    ur = model.eval_point([0, calibration_image_shape[1]], 'f').flatten()
    ll = model.eval_point([calibration_image_shape[0], 0], 'f').flatten()
    lr = model.eval_point([calibration_image_shape[0], calibration_image_shape[1]], 'f').flatten()

    leftbound   = int( np.ceil ( max(ul[1], ll[1])) )
    rightbound  = int( np.floor( min(ur[1], lr[1])) )
    topbound    = int( np.ceil ( max(ul[0], ur[0])) )
    bottombound = int( np.floor( min(ll[0], lr[0])) )

    # Don't crop left of lenslets (0, y) or above (x, 0)
    leftbound = max(leftbound, 0)
    topbound = max(topbound, 0)

    print leftbound, rightbound
    print topbound, bottombound

    return residuals[topbound:bottombound, leftbound:rightbound]


# ------------------------------------------------------------------
#                   DEBUGGING INFRASTRUCTURE
# ------------------------------------------------------------------

# Debug Image Files
FILE_GEOM_WHITE                = '00-GEOM-white.tif'
FILE_GEOM_WHITE_RECTIFIED      = '00-GEOM-white-rectified.tif'
FILE_GEOM_LOCALMAXIMA          = '00-GEOM-local-maxima.tif'
FILE_GEOM_CENTERS              = '00-GEOM-centers.tif'
FILE_CONVOLVED                 = '00-GEOM-convolved.tif'
FILE_TEMP                      = '00-GEOM-temp.tif'

# Misc. Files
FILE_LENSLETS_CACHE            = 'lenslet_grid_cache.npy'

def lenslet_debug_image(input_im, centers, lenslet_size, draw_grid = False):

    # Convert to RGB 
    out_im = np.tile(input_im,(3,1,1)).transpose(1,2,0)

    # (Optionally) draw the borders
    if (draw_grid):
        for i in range(0, out_im.shape[0], lenslet_size):
            out_im[i,:,:] = np.tile(np.array([0,0,1]), (out_im.shape[0],1))
        for i in range(0, out_im.shape[1], lenslet_size):
            out_im[:,i,:] = np.tile(np.array([0,0,1]), (out_im.shape[1],1))

        for i in range(lenslet_size, out_im.shape[0]+1, lenslet_size):
            out_im[i-1,:,:] = np.tile(np.array([0,0,1]), (out_im.shape[0],1))
        for i in range(lenslet_size, out_im.shape[1]+1, lenslet_size):
            out_im[:,i-1,:] = np.tile(np.array([0,0,1]), (out_im.shape[1],1))

    centers = np.cast['int16'](centers[:,0:2])

    # Remove centers that aren't in the image.
    centers = centers[np.nonzero(centers[:,0]>= 0)[0]]
    centers = centers[np.nonzero(centers[:,1]>= 0)[0]]
    centers = centers[np.nonzero(centers[:,0]< out_im.shape[0])[0]]
    centers = centers[np.nonzero(centers[:,1]< out_im.shape[1])[0]]

    # Draw centers as red.  (Swap indices here back no normal number col, row convention)
    out_im[centers[:,1], centers[:,0], :] = np.array([1,0,0])

    return out_im

# ------------------------------------------------------------------
#                   LIGHT FIELD CALIBRATION CODE
# ------------------------------------------------------------------

class CalibrationAlignmentMethods(object):

    CALIBRATION_ISOMETRY_ALIGNMENT = 1
    CALIBRATION_AFFINE_ALIGNMENT = 2
    CALIBRATION_CUBIC_ALIGNMENT = 3



def geometric_calibration(calibration_image, pixel_size, lenslet_pitch, alignment_method,
                          chief_rays = False, crop_center_lenslets = False):
    '''
    This geometric calibration routine takes either a radiometry or a
    chief ray light field.

    In the case where the lenslet circles show significant vignetting
    and appear more like "cat eyes", a more accurate calibration can
    be obtained by estimating the lenslet centers with a chief ray
    image rather than convolving the the radiometry image with a circular kernel.

    Paramaters:

      calibration_image - raw light field image (either radiometry or chief ray frame)

      pixel_size - size of a camera pixel in microns (taking relay lens maginifcation into account)

      lenslet_pitch - diameter of a lenslet in microns

      alignment_method - one of the CalibrationAlignmentMethods above

      chief_rays - set to 'True' if this is a chief ray image rather than a radiometry image.

      crop_center_lenslets - set to 'True' to limit the calibration fit to the
                             center lenslets in the array.  (Use in cases of severe vignetting.)

    Returns:

      an affine calibration matrix

    '''

    # Constants that tune algorithm performance
    ADJACENCY_TOLERANCE = 3 # Tolerance when fitting grid to lenslets.  3 works well in most situations.
    DEBUG = False           # Turn on debugging code.

    # -----------------
    # FIND LOCAL MAXIMA
    #
    approx_ppl = int(np.round(float(lenslet_pitch) / pixel_size))
    if approx_ppl < 1:
        raise RuntimeError("Error: arguments for lenslet pitch " + str(lenslet_pitch) + " and pixel size " + str(pixel_size) + " suggest that there is less than 1 pixel per lenslet!")

    from scipy.ndimage import filters
    if chief_rays == True:
        print '\t--> Calibrating lenslet centers using chief ray image image'

        # The chief ray image already has well isolated lenslet image
        # centers, and can be used directly in the code below.
        lf_tophat = calibration_image
        
    else:
        print '\t--> Calibrating lenslet centers using radiometry image'

        # Detect the local maxima using a difference of gaussians (to
        # detect lenslet-sized objects in the image) and a morphological
        # non-maximum suppression filter (to locate the center of each
        # lenslet).
        kernel = np.ones((approx_ppl, approx_ppl))
        for i in range(approx_ppl):
            for j in range(approx_ppl):
                radius = np.sqrt((i - approx_ppl/2.0)**2 + (j-approx_ppl/2.0)**2)
                if (radius >= approx_ppl/2.0-1):
                    kernel[i,j] = 0
        kernel = kernel / kernel.sum()  # Important: normalize the kernel

        # Convolve the top hat kernel with the image
        lf_tophat = filters.convolve(calibration_image, kernel)
        
    # Detect the local maxima using a difference of gaussians (to
    # detect lenslet-sized objects in the image) and a morphological
    # non-maximum suppression filter (to locate the center of each
    # lenslet).
    lf_nonmaxsup = nonmaxsup(lf_tophat, 10, 0.3)  # Threshold of 0.3    
    local_maxima = np.nonzero(lf_nonmaxsup)

    # Note: we use a [col, row] ordering in putative centers, which is
    # DIFFERENT than the natural [row, col] ordering commonly used to
    # index numpy matrices.  This is because we are dealing with image
    # data here, and we'll go crazy trying to keep track of when to
    # switch the indices.  Better to stay consistent.
    putative_centers = np.array([local_maxima[1],
                                 local_maxima[0]]).transpose()

    # Remove centers that are near the edge of the image (within 10
    # pixels or so).  It will be difficult to find good subpixel
    # estimates of the center of these spots.
    idxs = np.nonzero(putative_centers[:,0] > approx_ppl)
    putative_centers = putative_centers[idxs[0], :]
    idxs = np.nonzero(putative_centers[:,0] < calibration_image.shape[1] - approx_ppl)
    putative_centers = putative_centers[idxs[0], :]
    idxs = np.nonzero(putative_centers[:,1] > approx_ppl)
    putative_centers = putative_centers[idxs[0], :]
    idxs = np.nonzero(putative_centers[:,1] < calibration_image.shape[0] - approx_ppl)
    putative_centers = putative_centers[idxs[0], :]

    # We can optionally trim all lenslets outside the center region
    # (where no aperture vignetting occurs). Trimming the lenslets
    # outside of the 25% central region may significantly affect
    # calibration quality for high NA (~1.3) objectives.
    if crop_center_lenslets:
        new_centers = np.zeros(putative_centers.shape)
        count = 0;    
        for i in range(putative_centers.shape[0]):
            # use a circular central region of the lenlets for affine/cubic warp calibration
            print '--> Cropping only center lenslets for calibration'
            crop_percentage = 0.75
            r_max = calibration_image.shape[1]/2*crop_percentage
            x0 = calibration_image.shape[1]/2
            y0 = calibration_image.shape[0]/2
            x = putative_centers[i,0]
            y = putative_centers[i,1]
            r = np.sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0))
            if r < r_max:
                new_centers[count,:] = putative_centers[i,:]
                count += 1
                    
        putative_centers = new_centers[0:count,:]

    print "\t--> Found %d chief rays." % putative_centers.shape[0]

    # ---------------------------------
    # COMPUTE CONNECTED GRID OF CENTERS
    #
    print "\t--> Building 2D lenslet graph"

    # Refine peak location estimates with sub-pixel precision.
    subpixel_centers = refine_centers_subpixel(putative_centers, lf_tophat)

    # Build a grid of lenslets from the putative center estimates.
    lenslets = build_grid(subpixel_centers, approx_ppl, ADJACENCY_TOLERANCE)
    
    print '\t--> Created a lenslet grid.  Found %d lenslets.' % len(lenslets)

    # For debugging:
    # import scipy.io as io
    # mdict = dict()
    # mdict['chief_rays'] = lenslets[:,0:2]
    # mdict['lenslets'] = lenslets[:,2:]
    # if (chief_rays):
    #     io.savemat('/home/broxton/caldebug_chiefray.mat', mdict)
    # else:
    #     io.savemat('/home/broxton/caldebug_radiometry.mat', mdict)

    # --------------------------------------------------------------
    # SOLVE FOR THE MAPPING FROM CAMERA PIXEL SPACE TO LENSLET SPACE
    #
    # We are ready to solve for the best fit.  We'll use linear
    # least squares and a cubic or affine warp.
    #
    try:
        from lflib.calibration.models import AffineWarp, CubicWarp, IsometryWarp
        from lflib.calibration import LightFieldCalibration
        
        if alignment_method == CalibrationAlignmentMethods.CALIBRATION_AFFINE_ALIGNMENT:
            print '\t--> Computing Affine Warp'
            campixel_to_camlens = AffineWarp()
            
        elif alignment_method == CalibrationAlignmentMethods.CALIBRATION_CUBIC_ALIGNMENT:
            print '\t--> Computing Cubic Warp'
            campixel_to_camlens = CubicWarp()
            
        else:
            print '\t--> Computing Isometry Warp'
            campixel_to_camlens = IsometryWarp(lenslet_pitch, pixel_size)

        # Fit the data robustly using RANSAC
        bestdata = ransac(lenslets,
                          campixel_to_camlens,
                          40,                        # Minimum number of points needed for fit
                          20,                        # Number of iterations
                          0.5,                       # Threshold for inlier inclusion in the model
                                                     #     (0.5 == 1/2 a lenslet)
                          100)                       # Inliers required for 'good enough' model
        print '\t    RANSAC retained', bestdata.shape[0], '/', lenslets.shape[0], 'data points.'

        campixel_to_camlens.fit(bestdata)
        # campixel_to_camlens.save(os.path.join(output_filename))
        print campixel_to_camlens
        campixel_to_camlens.check_fit(bestdata)

        # Compute displacement vectors after the fit has been applied.
        # Measuring this deviation from the "ideal" model can provide
        # useful information about optical aberrations, such as the
        # defocus aberration introduced by the unknown distance
        # between the tube lens and the objective.
        projected_chiefrays = campixel_to_camlens.eval(lenslets, 'f')

        delta = lenslets[:, 2:] - projected_chiefrays
        calibration_residuals = np.zeros( (int(max(lenslets[:,3]) - 0.5)+1,
                                           int(max(lenslets[:,2]) - 0.5)+1),
                                          dtype = np.complex)
        for c in range(lenslets.shape[0]):
            calibration_residuals[int(lenslets[c,3]-0.5), int(lenslets[c,2]-0.5)] = delta[c,0] + delta[c,1] * 1j

        # DEBUGGING:  save out the residual images
        # print calibration_residuals.shape
        # distances = np.abs(calibration_residuals)
        # angles = np.angle(calibration_residuals)
        # save_image("_raw_distance.tif", distances.astype(np.float32))
        # save_image("_raw_angles.tif", angles.astype(np.float32))

        # Return the result
        return (campixel_to_camlens, calibration_residuals)

    except ValueError,e:
        raise RuntimeError('Calibration failed. RANSAC did not converge when fitting camera pixels to the lenslet grid.  Check your calibration file to make sure it has sufficient SNR & contrast.  The reported error was: ' + str(e))
        

