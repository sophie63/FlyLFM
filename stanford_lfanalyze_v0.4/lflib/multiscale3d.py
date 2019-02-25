# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2013 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
import time
import scipy.signal as sig
from scipy.ndimage.filters import median_filter

# R imports  
from rpy2.robjects.packages import importr
waveslim = importr("waveslim")
wavethresh = importr("wavethresh")
import rpy2.robjects as ro

# --------- Define R functions -------- 

ro.r( """
    "data_analytic_thresh" = function( wavelet_coefs, alpha=0.05, 
                                       max.level=6, verbose=FALSE,
                                       return.thresh=FALSE)
    {
      thresh = da.thresh(wavelet_coefs, alpha, max.level, verbose, return.thresh)
      thresh
    }

    "dwt_3d" = function( X_flattened, wavelet, J, vol_shape, undecimated=TRUE )
    {
      X = array(X_flattened, dim = c(vol_shape[2], vol_shape[1], vol_shape[3])) # x,y,z
      if (undecimated==TRUE)
      {
        wt = modwt.3d(X, wavelet, J)
        vol_inds = 0
      } else {
        J = ceiling(log(max(vol_shape),base=2))
        padded = pad_array(X, c(2^J, 2^J, 2^J))
        wt = wd3D(padded$Xpad)
        vol_inds = padded$vol_inds
      }
      list(wt=wt,vol_inds=vol_inds)
    }

    "pad_array" = function( array, final_shape )
    {
      out = array(0,dim=final_shape)
      n_pad = matrix(0,length(final_shape))
      for(i in 1:length(final_shape))
      {
        n_pad[i] = floor((final_shape[i]-dim(array)[i])/2)
      }
      vol_inds = c(n_pad[1],(n_pad[1]+dim(array)[1]-1),n_pad[2],(n_pad[2]+dim(array)[2]-1),n_pad[3],(n_pad[3]+dim(array)[3]-1))
      out[ vol_inds[1]:vol_inds[2], vol_inds[3]:vol_inds[4], vol_inds[5]:vol_inds[6] ] = array
      list(Xpad=out,vol_inds=vol_inds)
    }

    "hybrid_sure_thresh" = function( wavelet_coefs,
                                     max.level=4, 
                                     verbose=FALSE,
                                     seed=0)
    {
      thresh = hybrid.thresh(wavelet_coefs, max.level, verbose, seed) 
      thresh
    }

    "manual_thresh" = function( wavelet_coefs, value,
                                max.level=4,
                                hard=TRUE)
    {
      thresh = manual.thresh(wc, value=value, max.level=max.level, hard=hard)
      thresh
    }

    "idwt_3d" = function( wc, vol_inds=NULL, undecimated=TRUE )
    {
      if (undecimated==TRUE)
      {
        iwt = imodwt.3d(wc)
      } else {
        iwt = wr3D(wc)
        iwt = iwt[vol_inds[1]:vol_inds[2], vol_inds[3]:vol_inds[4], vol_inds[5]:vol_inds[6]]
      }
      iwt
    }

    "output_undecimated_wavelet_coefs" = function(wavelet_coefs, J)
    {
      output_coefs = list()
      for(i in 1:J)
      {
        output_coefs[[i]] = wavelet_coefs[[(i-1)*7 + 1]]
        for(j in 2:7)
        {
          output_coefs[[i]] = output_coefs[[i]] + wavelet_coefs[[(i-1)*7 + j]]
        }
      }
      output_coefs
    }

    "wavelet_coef_avg" = function(wavelet_coefs, undecimated=TRUE)
    {
      coef_sum = 0
      n = 0
      if(undecimated==TRUE)
      {
        n_coefs = length(wavelet_coefs)
        for( i in 1:n_coefs )
        {
          n = n + prod(dim(wavelet_coefs[[i]]))
          coef_sum = coef_sum + sum(sum(sum(wavelet_coefs[[i]])))
        }
      } else {
        n_levels = nlevelsWT(wavelet_coefs)
        for(i in 1:n_levels)
        {
          block_names = names(accessD(wavelet_coefs, level=i-1))
          n_blocks = length(block_names)
          for(j in 1:n_blocks)
          {
            coefs = accessD(wavelet_coefs, level=i-1, block=block_names[j])
            n = n + prod(dim(coefs))
            coef_sum = coef_sum + sum(sum(sum(coefs)))
          }
        }          
      }
      coef_avg = coef_sum/n
      coef_avg
    }

    "modify_wavelet_coefs" = function(wavelet_coefs, scale_factor, wavelet_mod, undecimated=TRUE)
    {
      update_SS = 0
      max_coef = 0
      if(undecimated==TRUE)
      {
        n_coefs = length(wavelet_coefs)
        for( i in 1:n_coefs )
        {
          update = scale_factor * wavelet_mod[[i]]
          wavelet_coefs[[i]] = wavelet_coefs[[i]] + update
          update_SS = update_SS + sum(sum(sum(update^2)))
          max_coef = max( max(wavelet_coefs[[i]]), max_coef )
        }
      } else {
        n_levels = nlevelsWT(wavelet_coefs)
        for(i in 1:n_levels)
        {
          block_names = names(accessD(wavelet_coefs, level=i-1))
          n_blocks = length(block_names)
          for(j in 1:n_blocks)
          {
            update = scale_factor*accessD(wavelet_mod, level=i-1, block=block_names[j])
            new_coefs = accessD(wavelet_coefs, level=i-1, block=block_names[j]) + update
            new_subarray = list(a=new_coefs, lev=i-1, block=block_names[j])
            wavelet_coefs = putD(wavelet_coefs, v=new_subarray)
            update_SS = update_SS + sum(sum(sum(update^2)))
          }
        }
      }
      update_norm = sqrt(update_SS)
      list( wavelet_coefs, update_norm, max_coef )
    }

    "sure_thresh" = function( wavelet_coefs,
                              max.level=4, 
                              hard=TRUE)
    {
      thresh = sure.thresh(wavelet_coefs, max.level, hard) 
      thresh
    }

    "universal_undecimated_thresh" = function( wavelet_coefs,
                                               max.level=4, 
                                               hard=TRUE)
    {
      thresh = universal.thresh.modwt(wavelet_coefs, max.level, hard)
      thresh
    }

   "wavelet_thresh" = function( wavelet_coefs, thresh=0, undecimated=TRUE, suppress=c(0) )
   {
     if(undecimated==TRUE)
     {
       n_coefs = length(wavelet_coefs)
       n_levels = floor(n_coefs/7)
       sigma = median(abs(wavelet_coefs$HHH1 - median(wavelet_coefs$HHH1)))/0.6745 #MAD

       if( length(thresh) == 1 ) {thresh = rep(thresh,n_coefs)}

       for( t in 1:n_levels )
       {
         for(j in 1:7)
         {
           idx = which( wavelet_coefs[[ (t-1)*7 + j ]] < thresh[[t]] )
           if(length(idx) > 0)
           {  
             wavelet_coefs[[ (t-1)*7 +j ]][idx] = 0
           } 
         }
       }

       for( s in suppress )
       {
         for(j in 1:7)
         {
           wavelet_coefs[[s*7 + j]] = 0*wavelet_coefs[[s*7 + j]] 
         }
       }
      
     } else {
       n_levels = nlevelsWT(wavelet_coefs)
       wavelet_coefs = threshold(wavelet_coefs, policy="manual", value=0, levels=0:7, by.level=T)
     }
     wavelet_coefs
   }

  """)

# ---- Import functions from R ----

# for testing
R_array = ro.r['array']

# wavelet transform functions
discrete_wavelet_transform_3d = ro.r['dwt_3d']
inverse_discrete_wavelet_transform_3d = ro.r['idwt_3d']

# wavelet thresholding functions
data_analytic_thresh = ro.r['data_analytic_thresh']
hybrid_sure_thresh = ro.r['hybrid_sure_thresh']
universal_undecimated_thresh = ro.r['universal_undecimated_thresh']
sure_thresh = ro.r['sure_thresh']

# Modify wavelet coefs by adding a scaled version of another set 
# of wavelet coefficients. 
modify_wavelet_coefs = ro.r['modify_wavelet_coefs']

# A manual threshold function
wavelet_thresh = ro.r['wavelet_thresh']

# Get average of wavelet coefficients
wavelet_coef_avg = ro.r['wavelet_coef_avg']
test_thresh = ro.r['threshold']

# For outputting coefficients for visualization
output_undecimated_wavelet_coefs = ro.r['output_undecimated_wavelet_coefs']

# These lines are required before we can covert numpy arrays into R arrays.
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

# ---- Callables ----

def multiscale_transform_3D(vol, vol_shape, wavelet_type='la8',transform_type="undecimated_wavelet"):

    J = int(np.ceil(np.log2(np.min(vol_shape))) )

    if transform_type == "DWT":
        coefs, vol_inds = discrete_wavelet_transform_3d( vol, wavelet=wavelet_type, 
                                                     J=J, vol_shape=vol_shape, 
                                                     undecimated=False)
    elif transform_type == "undecimated_wavelet":
        coefs, vol_inds = discrete_wavelet_transform_3d( vol, wavelet=wavelet_type, 
                                                     J=J, vol_shape=vol_shape, 
                                                     undecimated=True)
    elif transform_type == "pyramidal_median":
        vol = vol.reshape(vol_shape)
        coefs, vol_inds = hybrid_pyramidal_median_transform( vol, J=J, window_size=3) 
    else:
        raise ValueError("Specified transform type is not implemented.")

    return coefs, vol_inds

def inverse_multiscale_transform_3D( coefs, vol_inds, transform_type="undecimated_wavelet" ):
    if transform_type == "DWT":
        vol = inverse_discrete_wavelet_transform_3d(coefs,vol_inds,undecimated=False)
    elif transform_type == "undecimated_wavelet":
        vol = inverse_discrete_wavelet_transform_3d(coefs,vol_inds,undecimated=True)
    elif transform_type == "pyramidal_median":
        vol = inverse_hybrid_pyramidal_median_transform( coefs, vol_inds )
    else:
        raise ValueError("Specified transform type is not implemented.")

    return vol

def multiscale_coefficient_mean(coefs, transform_type="undecimated_wavelet"):
    if transform_type == "DWT":
        coef_mean = np.asarray(wavelet_coef_avg(coefs,'FALSE'))
    elif transform_type == "undecimated_wavelet":
        coef_mean = np.asarray(wavelet_coef_avg(coefs))
    elif transform_type == "pyramidal_median":
        coef_mean = pyramidal_median_coef_avg(coefs)
    else:
        raise ValueError("Specified transform type is not implemented.")

    return coef_mean

def multiscale_coefficient_update(coefs, update, update_rate, transform_type="undecimated_wavelet"):
    if transform_type == "DWT":
        new_coefs, update_norm, max_coef = modify_wavelet_coefs(coefs, update_rate, update, 'FALSE')
        update_norm = np.asarray(update_norm)[0]
    elif transform_type == "undecimated_wavelet":
        new_coefs, update_norm, max_coef = modify_wavelet_coefs(coefs, update_rate, update)
        update_norm = np.asarray(update_norm)[0]
    elif transform_type == "pyramidal_median":
        new_coefs, update_norm = modify_pyramidal_median_coefs(coefs, update, update_rate)
    else:
        raise ValueError("Specified transform_type is incorrect.")

    return new_coefs, update_norm

def multiscale_threshold( coefs, threshold=0.0, transform_type="undecimated_wavelet", 
                          suppress_scales=[0]):
    if transform_type == "DWT":
        thresholded_coefs = wavelet_thresh(coefs, threshold, 'FALSE')
    elif transform_type == "undecimated_wavelet":
        thresholded_coefs = wavelet_thresh(coefs, threshold, 'TRUE', suppress=suppress_scales)
    elif transform_type == "pyramidal_median":
        thresholded_coefs = threshold_pyramidal_median_coefs(coefs, threshold=threshold)
    else:
        raise ValueError("Specified transform type is not valid.")

    return thresholded_coefs

def output_multiscale_coefs(coefs, J, transform_type="undecimated_wavelet"):
    if transform_type == "DWT":
        raise NotImplementedError("DWT coefficient output not yet implemented.")
    elif transform_type == "undecimated_wavelet":
        output_coefs = output_undecimated_wavelet_coefs(coefs,J)
    elif transform_type == "pyramidal_wavelet":
        output_coefs = output_pyramidal_median_coefs(coefs,J)
    else:
        raise ValueError("Specified transform type is not valid.")

    return output_coefs

#-------- Variance stabilizing transforms -----------

def anscombe(x):
    return 2.0*np.sqrt(x + 3.0/8.0)

def inverse_anscombe(z):
    return (z/2.0)**2 - 3.0/8.0

def generalized_anscombe(x,mu,sigma,gain=1.0):
    return (2.0/gain)*np.sqrt(gain*x + (gain**2)*3.0/8.0 + sigma**2 - gain*mu)

def inverse_generalized_anscombe(z,mu,sigma,gain=1.0):
    return (1.0/gain)*(gain*y/2.0)**2 - gain*3.0/8.0 - (sigma**2)/gain + mu

# ------------ Multiscale Median Transforms ---------------

def multiscale_median_transform(x, J=4, initial_size=3):

    s = initial_size
    c = []; w = []
    c.append(x)
    for j in xrange(J):
        print "MMT level", j
        c.append( median_filter(c[j],size=s) )
        w.append( c[j] - c[j+1] )
        s *= 2

    return w.append( c[J] )

def hybrid_pyramidal_median_transform(vol, J=6, switch=2, window_size=3, verbose=True):

    # Pad volume out to a power of 2
    vol, unpad_idx = power2_pad_3D(vol)

    # initialize containers
    c = []; w = []
    c.append(vol)

    # structuring element for median filter
    struc_el = np.zeros((3,3,3)).astype(bool)
    struc_el[1,1,:] = 1
    struc_el[:,1,1] = 1
    struc_el[1,:,1] = 1

    s = 3
    # undecimated scales
    for j in xrange(switch):
        print "MMT level", j
        c.append( median_filter(c[j],size=s) )
        w.append( c[j] - c[j+1] )
        s *= 2

    # loop over J subbands
    for j in xrange(switch, J-1):
        if verbose:
            print "\t--> Pyramidal Median Transform level", j
        c.append( decimate_2x_3D( median_filter(c[j],size=s) ))
        c_upsample = upsample_2x_3D( c[j+1] )
        w.append( c[j] - c_upsample )
    w.append( c[-1] )

    return w, unpad_idx

def inverse_hybrid_pyramidal_median_transform(coefs, unpad_idx, J=6, switch=2):
    c = coefs[-1] # initialize with c_J 
    for j in reversed(xrange(switch,J-1)):
        c_upsample = upsample_2x_3D( c )
        c = coefs[j] + c_upsample 
    for j in reversed(xrange(switch)):
        c += coefs[j]
    c_unpadded = power2_unpad_3D( c, unpad_idx )
    return c_unpadded

def pyramidal_median_transform(vol, J=6, window_size=3, verbose=True):

    # Pad volume out to a power of 2
    vol, unpad_idx = power2_pad_3D(vol)

    # initialize containers
    c = []; w = []
    c.append(vol)

    # structuring element for median filter
    struc_el = np.zeros((3,3,3)).astype(bool)
    struc_el[1,1,:] = 1
    struc_el[:,1,1] = 1
    struc_el[1,:,1] = 1

    # loop over J subbands
    for j in xrange(J-1):
        if verbose:
            print "\t--> Pyramidal Median Transform level", j
        c.append( decimate_2x_3D( median_filter(c[j],footprint=struc_el) ))
        c_upsample = upsample_2x_3D( c[j+1] )
        w.append( c[j] - c_upsample )
    w.append( c[-1] )

    return w, unpad_idx

def inverse_pyramidal_median_transform(coefs, unpad_idx, J=6):
    c = coefs[-1] # initialize with c_J 
    for j in reversed(xrange(J-1)):
        c_upsample = upsample_2x_3D( c )
        c = coefs[j] + c_upsample 
    c_unpadded = power2_unpad_3D( c, unpad_idx )
    return c_unpadded

def pyramidal_median_coef_avg(coefs):
    n = 0
    coef_sum = 0.0
    for i in xrange(len(coefs)):
        subband = coefs[i]
        n += len(subband.flatten())
        coef_sum += np.sum(subband)
    return coef_sum/float(n)
    
def modify_pyramidal_median_coefs(coefs, update, update_rate=1.0):
    update_SS = 0.0
    for i in xrange(len(coefs)):
        update_subband = update_rate*update[i]
        update_SS += np.sum(update_subband**2)
        coefs[i] = coefs[i] + update_subband
    update_norm = np.sqrt(update_SS)
    
    return coefs, update_norm

def threshold_pyramidal_median_coefs(coefs, threshold):
    if type(threshold) == float:
        threshold = np.repeat(threshold, len(coefs))
    for i in xrange(len(coefs)):
        scale_coefs = coefs[i]
        scale_coefs[scale_coefs<threshold[i]] = 0.0
        coefs[i] = scale_coefs
    return coefs

def output_pyramidal_median_coefs(coefs, unpad_idx):
    vol_shape = coefs[0].shape
    for i in xrange(len(coefs)):
        while True:
            if coefs[i].shape[0] != vol_shape[0]:
                coefs[i] = upsample_2x_3D( coefs[i] )
                print vol_shape, coefs[i].shape
            else:
                break
    return coefs

# --------- Utility functions for resampling and padding 3D arrays ---------

def decimate_2x_3D( vol ):
    vol = sig.resample(vol, num=vol.shape[0]/2, axis=0)
    vol = sig.resample(vol, num=vol.shape[1]/2, axis=1)
    vol = sig.resample(vol, num=vol.shape[2]/2, axis=2)
    return vol

def upsample_2x_3D( vol ):
    vol = sig.resample(vol, num=vol.shape[0]*2, axis=0)
    vol = sig.resample(vol, num=vol.shape[1]*2, axis=1)
    vol = sig.resample(vol, num=vol.shape[2]*2, axis=2)
    return vol

def power2_pad_3D( vol ):
    n = 2**np.ceil(np.log2(vol.shape))
    padded_vol = np.zeros((n[0],n[1],n[2]))
    npad = [ np.floor((n[i] - vol.shape[i])/2) for i in range(3) ]
    idx = [ [npad[j], npad[j]+vol.shape[j]] for j in range(3) ]
    padded_vol[idx[0][0]:idx[0][1], idx[1][0]:idx[1][1], idx[2][0]:idx[2][1]] = vol
    return padded_vol, idx

def power2_unpad_3D( padded_vol, idx ):
    vol = padded_vol[idx[0][0]:idx[0][1], idx[1][0]:idx[1][1], idx[2][0]:idx[2][1]] 
    return vol

# ---- Testing ----

def test_wavelet_method_times(wavelet='la8', depth=None):
    """
    Function for testing speed of wavelet transforms 
    at different depths and wavelet functions.
    """
    # generate some data
    vol_shape = np.array([100,120,80])
    X = np.random.randn(vol_shape[0], vol_shape[1], vol_shape[2])
    X_flattened = np.array(X.ravel())

    # set parameters
    if depth is None:
        J = int(np.floor(np.log2(np.min(vol_shape))))
    else:
        J = depth
    print "J=",J
    J=4
    tic = time.time()
    wt, vol_inds = discrete_wavelet_transform_3d( X_flattened, wavelet, J, vol_shape, 'FALSE')
    toc = time.time() - tic
    print "Discrete wavelet transform took ", toc, " seconds."

    tic = time.time()
    iwt = inverse_discrete_wavelet_transform_3d(wt, vol_inds, 'FALSE')
    toc = time.time() - tic
    print "Discrete inverse wavelet transform took ", toc, " seconds."
    print "Shape of reconstructed volume:", np.asarray(iwt).shape

    return wt

def test_modify_wavelet_coefs( wavelet_coefs ):
    wavelet_coefs = test_wavelet_method_times()
    x, update_norm, max_coef = modify_wavelet_coefs( wavelet_coefs, 2, wavelet_coefs)
    print "norm:", np.asarray(update_norm)[0]

#-----------------------------------------------------------------------------
# Main. 

if __name__ == "__main__":

    tic = time.time()
    x = np.random.randn(157, 186, 60)
#    mmt = multiscale_median_transform(x)
    pmt, padded = hybrid_pyramidal_median_transform(x)
    x_recon = inverse_hybrid_pyramidal_median_transform(pmt,padded)

    print "MMT took", time.time() - tic, "seconds."
    1/0

    wt = test_wavelet_method_times()
    test_modify_wavelet_coefs(wt)

#EOF
