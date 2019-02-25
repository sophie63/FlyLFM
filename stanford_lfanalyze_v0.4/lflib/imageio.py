# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
import os
try:
    import cv2
except ImportError:
    import cv as cv2

def load_image(filename, dtype = None, normalize = False):
    """
    Load the image supplied by the user using OpenCV, and then
    immediately convert it to a numpy array.  The user may request
    that it be cast into a specific data type, or (in the case of
    floating point data) normalized to the range [0, 1.0].
    """
    if not os.path.exists(filename):
        raise IOError("File \"" + filename + "\" does not exist.")
    
    filetype = filename.split('.')[-1]
    if filetype.lower() == 'tif':
        from libtiff import TIFF
        tif = TIFF.open(filename, mode='r')

        # Each tiff directory contains one z slice!
        z_count = 0
        for zslice in tif.iter_images():

            # Handle Endian-ness conversion since pylibtiff doesn't do it automatically for some reason.
            if tif.IsByteSwapped():
                zslice = zslice.byteswap()
                
            shape = zslice.shape
            if z_count == 0:
                im = np.zeros((shape[0], shape[1], 1), dtype=zslice.dtype)
                im[:,:,0] = zslice
            else:
                im = np.concatenate((im, np.reshape(zslice, (shape[0], shape[1], 1))), axis=2)
            z_count += 1

        # If the tiff image has only one dimension, we squeeze it out of existence here.
        if z_count == 1:
            im = np.squeeze(im)

        del tif # Close the image
    elif filetype.lower() == 'jpg':
        # convert RGB to monochromatic
        im = np.asarray(cv2.imread(filename, cv2.CV_LOAD_IMAGE_GRAYSCALE))
    else:
        try:
            im = np.asarray(cv2.imread(filename, -1))
        except:
            im = np.asarray(cv2.LoadImage(filename, -1))
            im = np.asarray(im.ravel()[0][:]) # hack
            print "You are using an old version of openCV. Loading image using cv.LoadImage."

    if not im.shape:
        raise IOError("An error occurred while reading \"" + filename + "\"")

    # The user may specify that the data be returned as a specific
    # type.  Otherwise, it is returned in whatever format it was
    # stored in on disk.
    if dtype:
        im = im.astype(dtype)

    # Normalize the image if requested to do so.  This is only
    # supported for floating point images.
    if normalize :
        if (im.dtype == np.float32 or im.dtype == np.float64):
            return im / im.max()
        else:
            raise NotImplementedError
    else:
        return im


def save_image(filename, image, dtype = None):
    """
    Save the image to disk using OpenCV or libtiff.  The file format
    to use is inferred from the suffix of the filename.  OpenCV is
    used to write all file types except for tiff images.

    When saving tiffs, you may a 2D or 3D image into save_image().  A
    3D image will be saved as a tif stack automatically.
    """
    filetype = filename.split('.')[-1]

    # Test if there is no filetype
    if filename == filetype or len(filetype) > 3:
        raise IOError('Could not write file \'%s\'.  You must specify a file suffix.' % (filename))

    if os.path.dirname(filename) and not os.path.exists(os.path.dirname(filename)):
        raise IOError("Directory \"" + os.path.dirname(filename) +
                      "\" does not exist.  Could not save file \"" + filename + "\"")

    # If the user has specified a data type, we convert it here.
    if dtype != None:
        image = image.astype(dtype)

    # For now we transpose the data since it is stored in y,x,z order.
    # We can remove this later when we switch to z,y,x.
    if len(image.shape) == 3:
        image = np.transpose(image, (2,0,1))

    if filetype.lower() == 'tif':
        from libtiff import TIFF
        tif = TIFF.open(filename, mode='w')
        tif.write_image(image)
        del tif # flush data to disk
    elif len(image.shape) == 2:
        try:
            cv2.imwrite(filename, image)
        except:
            print "You are using an old version of openCV. Saving image using cv.SaveImage."
            try:
                cv2.SaveImage(filename, image)
            except:
                print 'There was an error saving the file', filename
    else:
        raise IOError('Saving 3-D image cubes to a \'%s\' file is not supported.  Please save your image as a \'tif\' stack.' % (filetype))
        

