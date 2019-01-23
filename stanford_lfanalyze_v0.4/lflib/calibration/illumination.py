import numpy as np
import numpy.linalg as linalg
import pyopencl as cl

from scipy.ndimage import filters
import scipy.io

from calibration.constants import *
from calibration.frames import AffineWarp
from calibration.io import *
from calibration.grid import buildGrid
from calibration.ransac import ransac
from calibration.graycode import  graycodeToBinary

# -------------------------------------------------------------------
# postProcessImage()
#
# This function cleans up the graycode image by (1) radiometrically
# calibrating it by dividing by the whitefield image, (2) 
#
# Returns a few different things:
#
# stmean - s x t image containing avg. brightness for each lenslet in the whitefield image.
# scatter - s*u x t*v image containing blurred scatter for this input image
# corrected - s*u x t*v image containing the radiometrically and geometrically corrected graycode image.
# 
def postProcessImage(graycode, white, uv_size, st_size):

    #cl_ctx = cl.create_some_context(False) # Allows you to select which device to use interactively
    devices = cl.get_platforms()[0].get_devices()
       
    # XXX: HACK!!  We select the first device here in a hard-coded
    # manner b/c I'm not sure how else to force OpenCL to use the GPU
    # instead of the CPU.  This should clearly be changed.
    cl_ctx = cl.Context(devices = [ devices[1] ])
    cl_queue = cl.CommandQueue(cl_ctx)
    mf = cl.mem_flags
    block_size = 16
    global_size = (int(np.ceil(float(st_size[0]) / block_size) * block_size),
                   int(np.ceil(float(st_size[1]) / block_size) * block_size))


    # STAGE 1 : COMPUTE STMEAN AND SCATTER IMAGES
    mod_postprocess_stage1 = cl.Program(cl_ctx, """
           __kernel void postprocess(unsigned short const usize, unsigned short const vsize,
                                     unsigned short const ssize, unsigned short const tsize,
                                     __global float* const graycode, __global float* const white,
                                     __global float* stmean, __global float* scatter) {

          // Calculate idx of the lenslet we are working on.
          unsigned int s = get_global_id(0);
          unsigned int t = get_global_id(1);
          unsigned idx = t * ssize + s;

          // We might be outside the reachable pixels. Don't do anything
          if( (s >= ssize) || (t >= tsize) )
            return;

          // Iterate over pixels in the input image.  Compute the mean.
          float mean_sum = 0;
          for (int v = 0; v < vsize; v++) {
            for (int u = 0; u < usize; u++) {
          
              // Sum up the pixels in the whitefield image, so that we
              // can later use this value as the mean for this disk.
              mean_sum += white[ ((t*vsize + v) * (usize*ssize)) + (s*usize + u) ];

              // Sum up the scatter for the graycode image.
              float scatter_sum = 0;
              for (int uu=0; uu<usize; uu++) {
                scatter_sum += graycode[ ((t*vsize + v) * (usize*ssize)) + (s*usize + uu) ];
              }
              for (int vv=0; vv<vsize; vv++) {
                scatter_sum += graycode[ ((t*vsize + vv) * (usize*ssize)) + (s*usize + u) ];
              }
              
              unsigned scatter_idx = ((t*vsize + v) * (usize*ssize)) + (s*usize + u);
              scatter[scatter_idx] = scatter_sum / (usize+vsize);
            }
          }


          // calculate offset into destination array
          stmean[idx] = mean_sum/(usize * vsize);
      }
      """).build()
     
    # Set up the buffers
    stmean = np.zeros(st_size, dtype="float32")
    scatter = np.zeros(graycode.shape, dtype="float32")
    graycode_buf = cl.Buffer(cl_ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=graycode)
    white_buf = cl.Buffer(cl_ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=white)
    stmean_buf = cl.Buffer(cl_ctx, mf.WRITE_ONLY, stmean.nbytes)
    scatter_buf = cl.Buffer(cl_ctx, mf.WRITE_ONLY, scatter.nbytes)

    # Execute the kernel on the GPU!!
    mod_postprocess_stage1.postprocess(cl_queue, global_size, (block_size, block_size),
                                       np.uint16(uv_size[1]), np.uint16(uv_size[0]),
                                       np.uint16(st_size[1]), np.uint16(st_size[0]),
                                       graycode_buf, white_buf, stmean_buf, scatter_buf)

    # Read back the result
    cl.enqueue_copy(cl_queue, stmean, stmean_buf).wait()
    cl.enqueue_copy(cl_queue, scatter, scatter_buf).wait()

    # Do a quick smoothing of the scatter image.  Marc's original code
    # actually smoothed across s,t space, whereas this filter operates
    # (mostly) in uv space.  I think the effect will probably be
    # similar, and it's much easier to call out to scipy real quick
    # here than to gin up another cuda kernel for the task.
    scatter = filters.uniform_filter(scatter, size=3, mode='constant')

    # ---

    # STAGE 2 : COMPUTE RADIOMETRICALLY CORRECTED GRAYCODE IMAGE AND GRADIENT IMAGES
    mod_postprocess_stage2 = cl.Program(cl_ctx, """
          __kernel void postprocess(unsigned short usize, unsigned short vsize,
                                    unsigned short ssize, unsigned short tsize,
                                    __global float* graycode, __global float* whitefield,
                                    __global float* blurred_scatter,
                                    float scatter_attenuation, float noise_floor,
                                    unsigned short margin, __global float* corrected) {
                                     
          // Calculate idx of the lenslet we are working on.
          unsigned int s = get_global_id(0);
          unsigned int t = get_global_id(1);

          // We might be outside the reachable lenslets. Don't do anything if this is the case.
          if( (s >= ssize) || (t >= tsize) )
              return;
              
          // Apply radiometric correction to each pixel based on
          // corresponding pixels in rectified fullwhite image.
          // Also subtract a fraction of mean intensity of disk,
          // to correct for pollution of nearby pixels, and
          // subtract from fullwhite to prevent lower corr_value
          // If fullwhite pixel is zero, zero data pixel too
          //
          for (int v = 0; v < usize; v++) {
            for (int u = 0; u < usize; u++) {
              // Pixel index [u,v,s,t]
              unsigned idx = ((t*vsize + v) * (usize*ssize)) + (s*usize + u);

              // Subtract the scatter value from the pixel
              float scatter = scatter_attenuation * blurred_scatter[idx];
              float value = graycode[idx] - scatter;

              // Zero out pixels below a noise floor. Also ensures
              // that value - scatter is non-negative.
              if (value < noise_floor) value = 0;

              float white = whitefield[idx];
              if (white < noise_floor) white = 0;

              // Zero a margin around microlenses, as above
              if (u < margin || u >= usize-margin || v < margin || v >= vsize-margin)
                white = 0;

              // Apply radiometric correction to the remaining pixels by dividing by the 
              // white value.  This should standardize the illumination across the image.
              float corr_value;
              if (white > 0) {
                corr_value = value / white;
                corrected[idx] = min(max(corr_value,0.0f),1.0f); // Clamp to the range [0, 1]
              } else {
                corrected[idx] = 0;
              }
            }
          }
      }
      """).build()
    
    # Set up the buffers
    corrected = np.zeros(graycode.shape, dtype="float32")
    corrected_buf = cl.Buffer(cl_ctx, mf.READ_WRITE, corrected.nbytes)

    # Execute the kernel on the GPU!!
    mod_postprocess_stage2.postprocess(cl_queue, global_size, (block_size, block_size),
                                       np.uint16(uv_size[1]), np.uint16(uv_size[0]),
                                       np.uint16(st_size[1]), np.uint16(st_size[0]),
                                       graycode_buf, white_buf, scatter_buf,
                                       np.float32(SCATTER_ATTENUATION), np.float32(NOISE_FLOOR),
                                       np.uint16(MARGIN), corrected_buf)

    # Read back the result
    cl.enqueue_copy(cl_queue, corrected, corrected_buf).wait()

    # ---

    # Apply an auto-levels like adjustment that stretches the
    # non-saturated intensity values!  This is done by computing the
    # histogram, throwing out the black and white saturated bins, and
    # normalizing the remaining bins.
    hist = np.histogram(corrected, range = [0, 1.0], bins=100)[0]
    hist = hist[1:-1]                              # Throw out saturated bins
    cdf = np.cumsum(hist / np.double(hist.sum()))  # Normalize, and compute CDF
    try:
        low = np.nonzero(cdf > 0.05)[0][0] / 100.0;
        high = np.nonzero(cdf < 0.95)[0][-1] / 100.0;
        print '\t        Normalizing to range: ', low, ', ', high
        corrected = (corrected - low) / (high-low)
        corrected[np.nonzero(corrected < 0.0)] = 0;
        corrected[np.nonzero(corrected > 1.0)] = 1.0;

    except IndexError:
        print ('WARNING: Insufficient image contrast found when trying to normalize',
               'the graycode image.  Skipping Normalization.')

    # ---

    # STAGE 3 : COMPUTE AVG. GRADIENTS
    mod_postprocess_stage3 = cl.Program(cl_ctx, """
          __kernel void postprocess(unsigned short usize, unsigned short vsize,
                                    unsigned short ssize, unsigned short tsize,
                                    float gradient_scale, unsigned short gradient_mincount,
                                    __global float* corrected, __global float* stcontrast) {

          // Calculate idx of the lenslet we are working on.
          unsigned int s = get_global_id(0);
          unsigned int t = get_global_id(1);

          // We might be outside the reachable lenslets. Don't do anything if this is the case.
          if( (s >= ssize) || (t >= tsize) )
              return;
              
          // Apply radiometric correction to each pixel based on
          // corresponding pixels in rectified fullwhite image.
          // Also subtract a fraction of mean intensity of disk,
          // to correct for pollution of nearby pixels, and
          // subtract from fullwhite to prevent lower corr_value
          // If fullwhite pixel is zero, zero data pixel too
          //
          float gradient_sum = 0;
          int gradient_count = 0;
          for (int v = 0; v < usize-1; v++) {
            for (int u = 0; u < usize-1; u++) {
              // Pixel index [u,v,s,t]
              unsigned idx = ((t*vsize + v) * (usize*ssize)) + (s*usize + u);
              unsigned idx_u1 = ((t*vsize + v) * (usize*ssize)) + (s*usize + u+1);
              unsigned idx_v1 = ((t*vsize + v+1) * (usize*ssize)) + (s*usize + u);

              // COMPUTE GRADIENTS
              float value = corrected[idx];

              // Ignore (zeroed) pixels (from radiometric calibration)
              if (value == 0) continue;
              float valueu1 = corrected[idx_u1];

              if (valueu1 > 0) {
                gradient_sum += fabs(valueu1-value);
                gradient_count++;
              }

              float valuev1 = corrected[idx_v1];
              if (valuev1 > 0) {
                gradient_sum += fabs(valuev1-value);
                gradient_count++;
              }
            }
          }

          // Ignore microlenses with few valid pixels (at edges)
          // Doesn't work perfectly !!
          unsigned idx = t * ssize + s;
          if (gradient_count > gradient_mincount)
            stcontrast[idx] = gradient_sum * gradient_scale / gradient_count;
          else
            stcontrast[idx] = 0;
      }
    """).build()

    # Set up the buffers
    stcontrast = np.zeros(st_size, dtype="float32")
    stcontrast_buf = cl.Buffer(cl_ctx, mf.WRITE_ONLY, stcontrast.nbytes)

    # Execute the kernel on the GPU!!
    mod_postprocess_stage3.postprocess(cl_queue, global_size, (block_size, block_size),
                                       np.uint16(uv_size[1]), np.uint16(uv_size[0]),
                                       np.uint16(st_size[1]), np.uint16(st_size[0]),
                                       np.float32(GRADSCALE), np.uint16(MINCOUNT),
                                       corrected_buf, stcontrast_buf)

    # Read back the result
    cl.enqueue_copy(cl_queue, stcontrast, stcontrast_buf).wait()

    # Return all the results!
    return (stmean, corrected, stcontrast)


# -------------------------------------------------------------------

def doPostProcess():

    # Load the lenslet warp for imaging
    campixel_to_camlens = AffineWarp()
    campixel_to_camlens.load(CALIBRATION_CAMPIXEL_TO_CAMLENS)
    
    # Load the eip-illumination whitefield image, which we will
    # use for calibration.
    print "\t--> Rectifying epi-illumination white calibration image."
    white = load_image(FILE_GEOM_WHITE)
    white_rectified = campixel_to_camlens.warpImage(white, DESIRED_LENSLET_SIZE, 'r', RECTIFIED_IMAGE_SIZE)

    # Compute some useful image sizes
    uv_size = (DESIRED_LENSLET_SIZE, DESIRED_LENSLET_SIZE)
    st_size = (white_rectified.shape[1]/DESIRED_LENSLET_SIZE,
               white_rectified.shape[0]/DESIRED_LENSLET_SIZE)

    print "\t    UV Size: ", uv_size
    print "\t    ST Size: ", st_size

    print '\t--> Radiometrically correcting graycode images.'
    for direction in ['h', 'v']:
        for bit in range(GRAYCODE_BITS):
            in_filename = (FILE_STRIPES_GRAYCODE + '-' + direction + '-' + str(bit) + '.png')
            out_filename = (FILE_POSTPROCESS_CORRECTED + '-' + direction + '-' + str(bit) + '.png')
            print '\t    ' + in_filename 

            # First we geomtrically rectify the graycode image
            graycode_im = load_image(in_filename)
            graycode_rectified = campixel_to_camlens.warpImage(graycode_im, DESIRED_LENSLET_SIZE,
                                                               'r', RECTIFIED_IMAGE_SIZE)

            # Then we send it off for post-processing and radiometric
            # calbration.
            (stmean, corrected, stcontrast) = postProcessImage(graycode_rectified, white_rectified,
                                                               uv_size, st_size)
            save_image(out_filename, corrected)

            # We only need to save out the stmean and stcontrast
            # images once for the "best" bit used when determining
            # where we have bad lenlet overlaps between the
            # illumination and imaging arrays.
            if (bit == 1):
                save_image(FILE_POSTPROCESS_STMEAN, stmean)
                save_image(FILE_POSTPROCESS_STCONTRAST + '-' + direction + '.png', stcontrast)

def doDecode():

    # Load stmean, and the "bestbit" stcontrast image, which is bit 1
    stmean = load_image(FILE_POSTPROCESS_STMEAN)
    stcontrast_s = load_image(FILE_POSTPROCESS_STCONTRAST + '-v.png')
    stcontrast_t = load_image(FILE_POSTPROCESS_STCONTRAST + '-h.png')

    # Create the 'stmask' image, which zeros out any lenslets with a
    # low gradient in bit 1 of the gray code.  These are camera
    # lenslets that probably fall on the the border between
    # illumination lenslets, and should therefore be ignored when
    # reading the gray code values.
    stmask = np.cast[np.float32]((stmean > MEANTHRESH) &
                                 (stcontrast_s > GRADTHRESH) & (stcontrast_t > GRADTHRESH))
    print '\n\t=============== Camera lenslet to projector pixel ==================='
    print '\t--> Saving stmask image.'
    save_image(FILE_DECODE_STMASK, stmask)

    coords = {}
    for direction in ['h', 'v']:
        print '\t--> Decoding and calibrating imaging pixels in the ', direction, ' direction.'
        # Load up the calibrated graycode images for this image dimension.
        images = [];
        graycode_image = np.zeros(ST_IMAGE_SIZE,dtype="uint16")
        for bit in range(GRAYCODE_BITS):
            in_filename = (FILE_POSTPROCESS_CORRECTED + '-' + direction + '-' + str(bit) + '.png')
            out_filename = (FILE_DECODE_THRESHOLD + '-' + direction + '-' + str(bit) + '.png')
            raw_image = load_image(in_filename)
            # For debugging:
            save_image(out_filename, np.cast['float32'](raw_image > 0.5))
            
            # Select the center lenslet location (this is the pixel
            # corresponding to the principal axis of the lenslet,
            # which is the one we can expect to be the most clear) and
            # apply the stmask.
            base_position = np.round(DESIRED_LENSLET_SIZE/2)
            principal_axis_img = stmask * raw_image[base_position::DESIRED_LENSLET_SIZE,
                                                    base_position::DESIRED_LENSLET_SIZE]

            # Threshold the result, and add these bits onto the
            # graycode image using the logical or operator...
            principal_axis_img = np.cast[np.uint16](principal_axis_img > 0.5)

            # For debugging
            #pa_filename = (FILE_DECODE_PAXIS + '-' + direction + '-' + str(bit) + '.png')
            #save_image(pa_filename, principal_axis_img)

            graycode_image |= (principal_axis_img << bit)

        coords[direction] = graycodeToBinary(graycode_image, GRAYCODE_BITS)
    out_filename = (FILE_DECODE_COORD)
    scipy.io.savemat(out_filename, {'graycode_coords_h':coords['h'],
                                    'graycode_coords_v':coords['v']})

    # Parse out the correspondences
    h = coords['h']
    v = coords['v']
    correspondences = np.zeros((h.shape[0]*h.shape[1],4))
    c = 0
    for s in range(h.shape[1]):
        for t in range(h.shape[0]):
            if (coords['v'][t,s] != 0 and coords['h'][t,s] != 0):
                correspondences[c,:] = np.array( [ float(s), float(t),
                                                 coords['v'][t,s], coords['h'][t,s]])
                c += 1

    # Get rid of extra zero entries in the matrix
    correspondences = correspondences[1:c,:]
    print '\t--> Found', c, 'correspondences.'

    # SOLVE FOR THE MAPPING FROM PIXEL SPACE TO LENSLET SPACE
    #
    # We are ready to solve for the best fit.  We'll use linear
    # least squares and a cubic warp.
    #
    try:
        camlens_to_projpixel = AffineWarp()
        bestdata = ransac(correspondences,
                          camlens_to_projpixel,
                          10,                           # Minimum number of points needed for fit
                          200,                          # Number of iterations
                          ILLUMINATION_LENSLET_SIZE/3,  # Threshold for inlier inclusion in the model
                          50)                           # Inliers required for 'good enough' model
        print '\t--> RANSAC retained',bestdata.shape[0],'/', correspondences.shape[0],'data points.'
        camlens_to_projpixel.fit(bestdata)
        camlens_to_projpixel.save(CALIBRATION_CAMLENS_TO_PROJPIXEL)
        print camlens_to_projpixel
        camlens_to_projpixel.checkFit(bestdata)
        
        # Save the projected camera pixels so that we can find the
        # projector pixel to projector lens transform in the next step.
        projected_campixels = bestdata[:,2:4]
        np.save(FILE_DECODE_PROJCENTERS, projected_campixels)
    except ValueError:
        print 'RANSAC did not converge when fitting imaging lenslets to projector pixels.'
    
def doCalibrateProjector():

    # Load up the positions of camera centers as determined by the
    # graycode decoding step.
    projected_campixels = np.load(FILE_DECODE_PROJCENTERS)

    # Camera microlenses are smaller than illumination lenslets; thus,
    # multiple microlenses will map to the same lenslet center.  These
    # nearly identical dots will confuse the cubic warp finder.  To
    # eliminate them, for all illuminated pixels, search for nearby
    # pixels.  If found, erase both and turn on their average
    # position.
    SEARCH_RADIUS = ILLUMINATION_LENSLET_SIZE / 2.0

    # First pass: combine nearby centers by averaging their locations.
    putative_centers = np.zeros(projected_campixels.shape)
    count = 0
    for i in range(projected_campixels.shape[0]):
        t = np.tile(projected_campixels[i,:], (projected_campixels.shape[0], 1))
        z = np.sqrt(np.sum((projected_campixels - t)**2,1))
        matches = projected_campixels[np.nonzero(z < SEARCH_RADIUS)[0], :]
        if matches.shape[0] > 1:
            count += 1

        centroid = np.sum(matches, 0)/matches.shape[0]
        putative_centers[i,:] = centroid
    print '\n\t=============== Projector pixel to projector lenslet ==================='
    print '\t--> Merged ', count, '/', putative_centers.shape[0],'centroids.'
    
    # Use the grid finding algorithm to line these lenslets up to a grid!
    lenslets = buildGrid(putative_centers,
                         ILLUMINATION_LENSLET_SIZE,
                         ILLUMINATION_ADJACENCY_TOLERANCE)
    print '\t--> Fit a grid to ', len(lenslets), '/', putative_centers.shape[0], ' putative matches.'

    # Save out the local maxima image for the projector (for debugging)
    illumination_maxima = np.zeros((np.round(projected_campixels[:,1].max()+10),
                                    np.round(projected_campixels[:,0].max()+10)),
                                   dtype = 'float32')
    for l in lenslets:
        illumination_maxima[l[1], l[0]] = 1.0
    save_image(FILE_PROJECTOR_LOCALMAXIMA, illumination_maxima)

    # SOLVE FOR THE MAPPING FROM PROJECTOR PIXEL SPACE TO PROJECTOR
    # LENSLET SPACE
    #
    # We are ready to solve for the best fit.  We'll use linear
    # least squares and a cubic warp.
    #
    projpixel_to_projlens = AffineWarp()
    bestdata = ransac(lenslets,
                      projpixel_to_projlens,
                      40,                        # Minimum number of points needed for fit
                      20,                        # Number of iterations
                      2.0/ILLUMINATION_LENSLET_SIZE,    # Threshold for inlier inclusion in the model
                      100)                       # Inliers required for 'good enough' model
    print '\t--> RANSAC retained',bestdata.shape[0],'/', lenslets.shape[0],'data points.'
    
    projpixel_to_projlens.fit(bestdata)
    projpixel_to_projlens.save(CALIBRATION_PROJPIXEL_TO_PROJLENS)
    print '\t--> Found the following cubic warp:'
    print projpixel_to_projlens
    projpixel_to_projlens.checkFit(bestdata)

