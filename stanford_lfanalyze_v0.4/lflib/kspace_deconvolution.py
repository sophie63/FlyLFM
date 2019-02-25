# import major libraries
import os.path
import numpy as np
import array
import struct
import math
import time

# import local modules

# Useful functions                               
pjoin = os.path.join

def lfdeconvolve_kspace( nt,                        # int
                         ns,                        # int
                         num_slices,                # int
                         spacingZ,                  # float not used yet
                         nu,                        # int
                         nv,                        # int
                         magnification,             # float
                         NA,                        # float
                         apertureimage_lf,          # array
                         centerZ      = 0.0,        # float
                         mediumRI     = 1.3333333,  # float
                         umPerLenslet = 125.0,        # float
                         focalstack = False ):      
    """
    Reconstruct a deconvolved volume stack from a single camera frame
    and save as a numpy array.
    """
    print '========================= K-SPACE DECONVOLUTION ================================='
    if 1==1:
        print '  Dimensions (nt, ns) ', nt, ns
        print '  num_slices  ', num_slices
        print '  nu, nv ', nu, nv
        print 'magnification', magnification
        print '          NA ', NA
        print '    mediumRI ', mediumRI
        print 'umPerLenslet ', umPerLenslet
        print '     centerZ ', centerZ
        print '    spacingZ ', spacingZ
    
    kz_fact = 3 # kz oversampling factor, choose an integer >=2
    
    # use different naming convention to match MATLAB code  
    numZ = num_slices
    numZk = kz_fact * numZ
    
    # generate synthetic light field data of an idea "delta function volume",
    # a volume with a single non-zero voxel at the center
    deltafunction_pinhole = np.array(np.zeros((nt,ns)),complex)
    deltafunction_pinhole[int(np.ceil(nt/2.0))][int(np.ceil(ns/2.0))] = 1 
    deltafunction_lf = np.tile(deltafunction_pinhole,(nv,nu))
        
    #Vk= grid_k_space( deltafunction_lf, #run this to see a reconstructed delta function volume
    Vk = grid_k_space( apertureimage_lf,    
                   nt,                  # int
                   ns,                  # int
                   numZ,                  # int
                   numZk,                 # int
                   spacingZ,              # float
                   nu,                    # int
                   nv,                    # int
                   magnification,         # float
                   NA,                    # float
                   mediumRI   = mediumRI, # float
                   umPerLenslet = umPerLenslet )    # float
    Vk_delta = grid_k_space( deltafunction_lf,
                   nt,                  # int
                   ns,                  # int
                   numZ,                  # int
                   numZk,                 # int
                   spacingZ,              # float
                   nu,                    # int
                   nv,                    # int
                   magnification,         # float
                   NA,                    # float
                   mediumRI   = mediumRI, # float
                   umPerLenslet = umPerLenslet )    # float
    Vk_linear = Vk.flatten(1);
    Vk_delta_linear = Vk_delta.flatten(1);

    if focalstack:
        print 'Creating focal stack only (no k-space deconvolution)'
    else:
        # deconvolve by dividing 3D k-space by the (non-zero! elements of)
        # k-space "PSF", which is the 3D k-space from a delta function LF.
        Vk_linear[np.abs(Vk_delta_linear)>1e-12] =  Vk_linear[np.abs(Vk_delta_linear)>1e-12]/abs(Vk_delta_linear[np.abs(Vk_delta_linear)>1e-12])
    
    # reshape linear-indexing back to 3D indexing 
    Vk = Vk_linear.reshape(Vk.shape,order='F')
    
    # 3D IFFT back to get deconvolved volume
    V = np.abs(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(Vk),Vk.shape)));
    
    # crop out extra z-slices in volume from kz-oversampling
    V = V[:,:,V.shape[2]/2-numZ/2+1-1:V.shape[2]/2+numZ/2];

    # calculate z-spacing
    u = nu-1; v = nv-1;
    theta, phi = uv_to_angle(u, v, nu, nv, NA, mediumRI)
    #theta=-1*np.arcsin((float(v)-(nv-1.0)/2.0)/((nv-1.0)/2.0)*NA/mediumRI);
    print "max theta=",theta/3.14*180, 'degrees'
    Tx = np.float(umPerLenslet) / magnification
    Tz = np.abs(np.float(Tx) / np.tan(theta)) / 2
    print 'k-space output z-spacing should be ', Tz ,'um per slice', ' x-spacing: ', Tx, 'um'

    return V.astype(np.float32)
    
    # reorder dimensions to make output consistent with previous MLSL code
    #V=np.transpose(V,(1,0,2))
    #volume_stack=V.flatten(1)
    #vol       = np.array( volume_stack )
    #vol_shape = ( num_slices, ns, nt )
    #return vol, vol_shape

#  Convert from discrete ray angle direction to an angle.  The
#  lenslet class implements this use the Abbe sine correction.
#  Returns a tuple (theta, phi) of angles relative to the optical
#  axis of the lenslet.
# Compute theta and phi based using the Abbe sin condition
# (see Marc's 2009 paper & the Oldenberg reference therein
# for details.) 
def uv_to_angle(u, v, nu, nv, NA, n):
    phi   = np.arcsin((float(v)-(nv-1.0)/2.0)/((nv-1.0)/2.0)*NA/n);
    theta = np.arcsin((float(u)-(nu-1.0)/2.0)/((nu-1.0)/2.0)*NA/n);
    return theta, phi
    
def grid_k_space( apertureimage_lf,
                   nt,                  # int
                   ns,                  # int
                   numZ,                  # int
                   numZk,                 # int
                   spacingZ,              # float
                   nu,                    # int
                   nv,                    # int
                   magnification,            # float
                   NA,                       # float
                   mediumRI     = 1.3333333, # float
                   umPerLenslet = 125. ):    # float
    #print "starting grid_k_space()"
    n = mediumRI
    
    # determine zmax from the max theta and phi. All z-indices will be scaled to zmax.
    u = nu-1; v = nv-1;
    theta, phi = uv_to_angle(u, v, nu, nv, NA, n)
    #print "max theta=", theta/3.14*180, "max phi=", phi/3.14*180
    x = float(ns)/2; y = float(nt)/2;
    zmax = (x*np.tan(theta) + y*np.tan(phi));
    zmax *= -1 # reverse z-direction
    
    # precompute meshgrids which are used in the loop.
    X, Y = np.meshgrid([float(x) for x in range(-int(np.floor(ns/2.0)),int(np.ceil(ns/2.0)))],
                        [float(x) for x in range(-int(np.floor(nt/2.0)),int(np.ceil(nt/2.0)))])
    xind, yind = np.meshgrid([int(x) for x in range(ns)],[int(y) for y in range(nt)])
    xind = xind.flatten(1); yind = yind.flatten(1) #used to convert zindex to linear array (:)
    
    # initialize 3D k-space with zeros
    Vk = np.zeros((nt,ns,numZk),complex)
    Vk_linear = Vk.flatten(1) 
    
    for v in range(nv):
        for u in range(nu):
            # Grab the current sub-aperture image, perform 2D FFT to get 2D slice to be
            # put in 3D k-space
            im_pinhole = np.array(apertureimage_lf[v*nt:v*nt + nt, u*ns:u*ns + ns],complex)
            #print 'debug', ns, nt, ns*nt, im_pinhole.shape, nu, nv, numZ, numZk
            k_slice = np.fft.fftshift( np.fft.fft2(np.fft.fftshift(im_pinhole) )) # FFT takes ~60% of loop execution time!
            k_slicelinear = k_slice.flatten(1); #linear-indexed version of Vk
            
            theta, phi = uv_to_angle(u, v, nu, nv, NA, n)
            
            # the 2D slice (k_slice) needs to be tilted by (theta,phi) before 
            # being put into 3D k-space (Vk_linear). compute the z-indices using
            # nearest-neighbor (in z only) resampling.
            z = (X*np.tan(theta) + Y*np.tan(phi))/zmax; #(-1,1)
                # (z+1)*(nzk-1)/2
            zindex = np.array(np.round((z+1.0)*(float(numZk)-1.0)/2.0),int); #-numZk/2:numZk/2+1
            zindexlinear = yind + nt*(xind) + ns*nt*(zindex.flatten(1));
            
            # actually (additively) insert 3D slice into 3D k-space
            Vk_linear[zindexlinear]+=k_slicelinear[:] #takes 20% of loop execution time
    
    # reshape linear-indexing to 3D indexing
    Vk=Vk_linear.reshape(Vk.shape,order='F')

    #print "ending grid_k_space()"   
    
    # return 3D k-space
    return Vk
#-----------------------------------------------------------------------------#  
#EOF

