# Frame of reference manager for scope.py
#
# We manage the following frame of reference here:
#
# - Camera pixels to camera lenslets
# - Camera lenslets to projector lenslets
# - Projector pixels to projector lenslets
# - Camera lenslets to x,y positions in micrometers in the native plane.
import numpy as np
from calibration.constants import *
import os, math, sys

# ------------------------------------------------------------------
#                             FRAME NAMES
# ------------------------------------------------------------------

CAMERA_PIXEL_FRAME = 1
PROJECTOR_PIXEL_FRAME = 2
CAMERA_LENSLET_FRAME = 3
PROJECTOR_LENSLET_FRAME = 4
SAMPLE_MICRON_FRAME = 5

# ------------------------------------------------------------------
#                      FRAME MANAGER CLASS
# ------------------------------------------------------------------

# FrameManager Class
#
# This class helps to keep track of the various frames of reference
# commonly used in the uScope application.
class FrameManager(object):

    def __init__(self):
        self.campixel_to_camlens = AffineWarp()
        self.projpixel_to_projlens = AffineWarp()
        self.camlens_to_projpixel = AffineWarp()

    # For the frame manager to reload the calibration files
    def reload_calibration(self, root_directory):
        try:
            self.campixel_to_camlens.load(os.path.join(root_directory, CALIBRATION_CAMPIXEL_TO_CAMLENS))
#            self.projpixel_to_projlens.load(os.path.join(root_directory, CALIBRATION_PROJPIXEL_TO_PROJLENS))
#            self.camlens_to_projpixel.load(os.path.join(root_directory, CALIBRATION_CAMLENS_TO_PROJPIXEL))
        except IOError:
            print 'ERROR: could not load calibration file!'
            sys.exit(1)
        self.camlens_to_sample = None

    # For the frame manager to reload the calibration files
    def save_calibration(self, directory):
        try: 
            self.campixel_to_camlens.save(os.path.join(directory,
                                                       os.path.basename(CALIBRATION_CAMPIXEL_TO_CAMLENS)))
            self.projpixel_to_projlens.save(os.path.join(directory,
                                                         os.path.basename(CALIBRATION_PROJPIXEL_TO_PROJLENS)))
            self.camlens_to_projpixel.save(os.path.join(directory,
                                                        os.path.basename(CALIBRATION_CAMLENS_TO_PROJPIXEL)))
        except IOError:
            print 'WARNING: could not load calibration files.  The system needs to be calibrated.'
        self.camlens_to_sample = None

    # Coordinate transform method
    def transform_coordinates(self, coords, src_frame, dst_frame):

        if (dst_frame == CAMERA_PIXEL_FRAME):
            if (src_frame == PROJECTOR_PIXEL_FRAME):
                return self.transform_coordinates(self.camlens_to_projpixel.reverse(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)
            if (src_frame == PROJECTOR_LENSLET_FRAME):
                return self.transform_coordinates(self.projpixel_to_projlens.reverse(coords),
                                                  PROJECTOR_PIXEL_FRAME, dst_frame)
            if (src_frame == CAMERA_LENSLET_FRAME):
                return self.campixel_to_camlens.reverse(coords)
            if (src_frame == SAMPLE_MICRON_FRAME):
                return self.transform_coordinates(self.camlens_to_sample.reverse(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)

        if (dst_frame == PROJECTOR_PIXEL_FRAME):
            if (src_frame == CAMERA_PIXEL_FRAME):
                return self.transform_coordinates(self.campixel_to_camlens.forward(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)
            if (src_frame == PROJECTOR_LENSLET_FRAME):
                return self.projpixel_to_projlens.reverse(coords)
            if (src_frame == CAMERA_LENSLET_FRAME):
                return self.camlens_to_projpixel.forward(coords)
            if (src_frame == SAMPLE_MICRON_FRAME):
                return self.transform_coordinates(self.camlens_to_sample.reverse(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)

        if (dst_frame == CAMERA_LENSLET_FRAME):
            if (src_frame == CAMERA_PIXEL_FRAME):
                return self.campixel_to_camlens.forward(coords)
            if (src_frame == PROJECTOR_LENSLET_FRAME):
                return self.transform_coordinates(self.projpixel_to_projlens.reverse(coords),
                                                  PROJECTOR_PIXEL_FRAME, dst_frame)
            if (src_frame == PROJECTOR_PIXEL_FRAME):
                return self.camlens_to_projpixel.reverse(coords)
            if (src_frame == SAMPLE_MICRON_FRAME):
                return self.camlens_to_sample.reverse(coords)

        if (dst_frame == PROJECTOR_LENSLET_FRAME):
            if (src_frame == CAMERA_PIXEL_FRAME):
                return self.transform_coordinates(self.campixel_to_camlens.forward(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)
            if (src_frame == CAMERA_LENSLET_FRAME):
                return self.transform_coordinates(self.camlens_to_projpixel.forward(coords),
                                                  PROJECTOR_PIXEL_FRAME, dst_frame)
            if (src_frame == PROJECTOR_PIXEL_FRAME):
                return self.projpixel_to_projlens.forward(coords)
            if (src_frame == SAMPLE_MICRON_FRAME):
                return self.transform_coordinates(self.camlens_to_sample.reverse(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)

        if (dst_frame == SAMPLE_MICRON_FRAME):
            if (src_frame == CAMERA_PIXEL_FRAME):
                return self.transform_coordinates(self.campixel_to_camlens.forward(coords),
                                             CAMERA_LENSLET_FRAME, dst_frame)
            if (src_frame == CAMERA_LENSLET_FRAME):
                return self.camlens_to_sample.forward(coords)
            if (src_frame == PROJECTOR_LENSLET_FRAME):
                return self.transform_coordinates(self.projpixel_to_projlens.reverse(coords),
                                                  PROJECTOR_PIXEL_FRAME, dst_frame)
            if (src_frame == PROJECTOR_PIXEL_FRAME):
                return self.transform_coordinates(self.camlens_to_projpixel.reverse(coords),
                                                  CAMERA_LENSLET_FRAME, dst_frame)

