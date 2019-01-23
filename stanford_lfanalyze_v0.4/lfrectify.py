#!/usr/bin/env python2.7

# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

# rectify.py
#
# Usage: rectify.py <lightfield_image.{tif,png,etc}> [--pixels-per-lenslet <ppl>]
#
# This script simply applies a rectification from a
# campixel_to_camlens.warp file to a single light field image.  The
# resulting rectified image should have the lenslets aligned with the
# horizontal and vertical dimensions of the image.  You can optionally
# specify the number of pixels per lenlet you would like in the output
# image, otherwise this value is computed for you based on the input
# imagery and the warp file.

from lflib.imageio import load_image, save_image
import sys, os, math
import numpy as np

if __name__ == "__main__":

    # Parse command line options
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-p", "--output-pixels-per-lenslet", dest="output_pixels_per_lenslet",
                      help="Specify the number of pixels per lenslet in the output image.")
    parser.add_option("-o", "--output-file", dest="output_filename",
                      help="Specify the output filename.")
    parser.add_option("-c", "--calibration-file", dest="calibration_file",
                      help="Specify the calibration file to use for rectification.")
    parser.add_option('-s', "--subaperture",
                      action="store_true", dest="subaperture", default=False,
                      help="Save out the light field image as tiled subapertures.")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        print 'You must supply at least one light field image to rectify.\n'
        parser.print_help()
        sys.exit(1)

    print 'Rectifying', len(args), 'images.'

    for filename in args:

        # Default calibration filename has a *.lfc suffix, and the same prefix
        if not options.calibration_file:
            fileName, fileExtension = os.path.splitext(filename)
            calibration_file = fileName + '.lfc'
        else:
            calibration_file = options.calibration_file

        # Default output filename has a -RECTIFIED suffix
        if not options.output_filename:
            fileName, fileExtension = os.path.splitext(filename)
            output_filename = fileName + '-RECTIFIED' + fileExtension
        else:
            output_filename = options.output_filename
        
        # Load the calibration data
        from lflib.calibration import LightFieldCalibration
        lfcal = LightFieldCalibration.load(calibration_file)

        # Rectify the image
        im = load_image(filename, normalize = False)
        input_dtype = im.dtype

        # Perform dark frame subtraction
        im = lfcal.subtract_dark_frame(im)

        # Rectify the image
        rectified_lf = lfcal.rectify_lf(im)


        # Optionally reformat the image so that sub-aperturs are tiled, rather than lenslet images.
        from lflib.lightfield import LightField
        if (options.subaperture):
            im = rectified_lf.asimage(LightField.TILED_SUBAPERTURE)
            print '\t--> Saving ', output_filename, 'as tiled sub-aperture image.'
        else:
            im = rectified_lf.asimage(LightField.TILED_LENSLET)
            print '\t--> Saving ', output_filename, 'as tiled lenslet image.'

        print im.max()
        print im.min()

        #im = im/im.max()*255
        #print im.max()
        #print im.min()
        #save_image(output_filename, im.astype(np.uint8), dtype=np.uint8)
        save_image(output_filename, im, dtype=np.uint16)
