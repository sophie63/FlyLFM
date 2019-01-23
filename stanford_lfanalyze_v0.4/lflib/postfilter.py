import numpy as np
from lflib.imageio import save_image

def remove_light_field_grid(vol, supersample_factor):
  
    # Constants -- control the size and threshold of the spectral median
    # filter.  You can probably leave these at their default values...
    patch_size = 20;
    threshold = 3;

    # Compute fourier spectra of each slice in the stack
    spectrum = np.zeros(vol.shape, dtype=np.complex)
    for i in range(vol.shape[2]):
        spectrum[:,:,i] = np.fft.fftshift(np.fft.fft2( vol[:,:,i] ))

    # Iterate over slices
    for slice_num in range( vol.shape[2] ):

        # Save the spetrum near the origin (DC + low frequency components).
        # We will replace this portion of the spectrum after all the
        # filtering below.
        #
        # Making a copy here is important because otherwise origin
        # simply points at the data in spectrum, rather than making a
        # deep copy.
        origin = spectrum[ (vol.shape[0]/2-patch_size):(vol.shape[0]/2+patch_size), 
                           (vol.shape[1]/2-patch_size):(vol.shape[1]/2+patch_size), 
                           slice_num ].copy();

        # STEP 1: Filter out the lattice of spots corresponding to the
        # lenslet spacing.  
        dx = vol.shape[1]/supersample_factor;
        dy = vol.shape[0]/supersample_factor;
        for x in np.arange(dx, vol.shape[1], dx):
            for y in np.arange(dy, vol.shape[0], dy):

                # Don't mess with the DC + low frequency components!
                #if x == vol.shape[1]/2 and y == vol.shape[0]/2:
                #    print 'SKIPPING'
                #    continue
                
                patch = spectrum[ np.round(y-patch_size/2):np.round(y+patch_size/2),
                                  np.round(x-patch_size/2):np.round(x+patch_size/2),
                                  slice_num ];
                md = np.median(np.abs(patch));
                patch[np.nonzero(np.abs(patch) / md > threshold)] = md;
                spectrum[ np.round(y-patch_size/2):np.round(y+patch_size/2), 
                          np.round(x-patch_size/2):np.round(x+patch_size/2),
                          slice_num ] = patch;

        # Restore DC + low frequency components that might have been filtered
        # out by the above operations
        spectrum[ (vol.shape[0]/2-patch_size):(vol.shape[0]/2+patch_size), 
                  (vol.shape[1]/2-patch_size):(vol.shape[1]/2+patch_size), 
                  slice_num ] = origin;

    # Inverse Fourier transform the cleaned up spectrum.  Return result!
    result = np.zeros_like(vol)
    for i in range(vol.shape[2]):
        result[:,:,i] = np.abs(np.fft.ifft2(np.fft.fftshift( spectrum[:,:,i] )));
        
    return result
  
  
