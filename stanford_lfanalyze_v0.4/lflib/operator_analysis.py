import os,sys
import h5py
import numpy as np

# lflib imports
from lflib.iterative_deconvolution import LightFieldOperator
from lflib.calibration import LightFieldCalibration
from lflib.volume import LightFieldProjection
from lflib.optics import compute_light_field_psf
from lflib.optics import LensletArray
from lflib.imageio import load_image, save_image

# PyCuda imports
import pycuda.autoinit 
import pycuda.gpuarray as gpuarray
import scikits.cuda.fft as cu_fft

# Scipy imports
from scipy.sparse.linalg import LinearOperator
from scipy.stats import ks_2samp

#-------------------------------------------------------------------------------------------------
# Classes and functions for evaluating resolution at the camera sensor
#-------------------------------------------------------------------------------------------------

def noisy_sensor_resolution(p1,p2,error_rate=30):
    """
    This function takes two points in 3D space, p1 and p2, projects 
    each onto the camera sensor (noiselessly), adds Poisson noise with 
    rate parameter error_rate to each, and runs a KS test to determine 
    if the two distributions are distinguishable.
    """
    # set lenslet array parameters
    nu = 27
    nv = 27
    ns = 21
    nt = 21
    ulens_pitch = 125
    ulens_focal_length = 2426
    objective_magnification = 20
    objective_na = 0.5
    medium_index = 1.33

    # Construct lenslet array object
    lenslet_array = LensletArray(nu, nv, ns, nt,
                                 ulens_pitch, ulens_focal_length, ulens_focal_length,
                                 objective_magnification, objective_na, medium_index,
                                 ulens_fill_factor = 1.0, pixel_fill_factor = 1.0,
                                 circular_ulens_profile = False, 
                                 center_wavelength = 509)  # Units: nanometers

    # Input list with (intensity,x,y,z,num_lenslets_in_psf,lenslet_array,wavelength_nm)
    # to compute_light_field_psf; wavelength currently fixed at 510nm with intensity = 1.0.
    psf0 = compute_light_field_psf( None, 1.0, p1[0], p1[1], p1[2], ns, lenslet_array, 510 )
    psf1 = compute_light_field_psf( None, 1.0, p2[0], p2[1], p2[2], ns, lenslet_array, 510 )

    # Add gaussian noise (making poisson intensity >30 assumption)
    # The shot noise variance for each nonzero pixel should be linearly 
    # related to the mean intensity of the corresponding pixel in p1.
    noise = psf0 * np.random.normal(loc=0.0, scale=1.0, size=psf0.shape)
    signal = psf1 - psf0 + np.random.normal(loc=0.0, scale=1.0, size=psf0.shape)

    # log likelihood ratio on continuous data (based on poisson shot noise)
    l0 = 2*psf0
    la = psf1 + psf0
    logL = np.sum( la*(np.log(la) - np.log(l0) - 1.0) + l0 )

    # log likelihood ratio on discrete (16-bit) data (based on poisson shot noise)
    psf_max = 2*np.max(psf0)
    psf0_discrete = (65535*(psf0/psf_max)).astype(np.uint16)
    psf1_discrete = (65535*(psf1/psf_max)).astype(np.uint16)
    l0 = 2.0*psf0_discrete
    la = psf1_discrete + psf0_discrete
    log_la = np.log(la); log_la[np.where(log_la==-np.inf)[0]]=0.0
    log_l0 = np.log(l0); log_l0[np.where(log_l0==-np.inf)[0]]=0.0
    logL_discrete = np.sum( la*(log_la - log_l0 - 1.0) + l0 )

    # save 16-bit pngs
    save_image('/home/logan/Documents/Results/Resolution/sensor/psf0.png',psf0_discrete)
    save_image('/home/logan/Documents/Results/Resolution/sensor/psf1.png',psf1_discrete)
    
    # KS test
    ks, pval = ks_2samp( signal.flatten(), noise.flatten() )
    print "KS statistic:",ks
    print "KS p-value:",pval
    print "log Likelihood ratio:",logL
    print "Discrete log Likelihood ratio:",logL_discrete
    return ks,pval,logL,logL_discrete

#-------------------------------------------------------------------------------------------------
# Classes and functions for evaluating resolution in volumes
#-------------------------------------------------------------------------------------------------

class CovarianceLinearOperator(object):
    """
    Applies A'A to a vector if cov_type=="geometric" or to a delta function in a volume if the 
    cov_type=="wave", where A is the optical system linear operator.

    A'A applied to a point in a volume results in the volumetric PSF for that point.
    """
    def __init__(self, A, cov_type="geometric", vol_shape=None):
        self.A = A
        self.cov_type = cov_type
        try:
            self.shape = np.prod(vol_shape)
            self.vol_shape = vol_shape
        except:
            print "For a geometric model, you must supply the shape of the volume the covariance matrix refers to as a tuple (y,x,z)."

    def matvec(self, x):
        if self.cov_type == "geometric":
            matvec = self.A.rmatvec(self.A.matvec(x)) 
        elif self.cov_type == "wave":
            forward = compute_light_field_psf(*x) # x is a list with (intensity,x,y,z,num_lenslets_in_psf,lenslet_array,wavelength_nm)
            matvec = self.A.rmatvec(forward) 
        return matvec

    def rmatvec(self, x):
        return self.matvec(x) # operator is symmetric

#-------------------------------------------------------------------------------------------------

def get_system_operator_pca(Cov,K=3):
    """ 
    Get the K largest eigenvalues and vectors of sparse linear operator A using Lanczos iteration (ARPACK).
    Note, this may be incredibly slow for large A. TODO: add iteration control. 
    """
    from scipy.sparse.linalg import eigsh
    return eigsh(Cov, K, which='LM')

def get_psf_vol(point, Cov, model_type="geometric", raydb=None, lenslet_array=None):
    """
    Get psf of a delta function at 'point' (y,x,z) in a volume by either:
    -- if model_type=="geometric":
         Apply the covariance operator Cov = A'A to the voxel indexed 
         by the voxel coordinates supplied in 'point'.
    -- if model_type=="wave":
         Use wave optics model to estimate sensor psf of a delta function in 
         the volume at 'point', then apply the appropriate A' to this (Ax). 
    """
    if model_type=="geometric":
        # generate volume with single voxel "on"
        vol = np.zeros(raydb.nvoxels, dtype=np.float32).reshape(raydb.ny, raydb.nx, raydb.nz)
        vol[point] = 1.0*raydb.supersample_factor**3
        vol_vec = vol.flatten()

        # get psf
        psf_vec = Cov.matvec(vol_vec)
        psf_vol = np.reshape(psf_vec, (raydb.ny, raydb.nx, raydb.nz))

    elif model_type=="wave":
        # get number of lenslets in aperture for splat
        objective_theta = np.arcsin(lenslet_array.na / lenslet_array.medium_index) 
        aperture_diameter = np.abs( 2*point[2]*np.tan(objective_theta) ) 
        num_lenslets_in_aperture = int(np.ceil(aperture_diameter / (lenslet_array.pitch / lenslet_array.magnification))) 

        # generate point list information for wave optics covariance operator 
        # (this contains arguments [intensity,x,y,z,num_lenslets_in_psf,lenslet_array,wavelength_nm])
        # and get psf
        point_list = [1.0, point[1], point[0], point[2], num_lenslets_in_aperture, lenslet_array, 510] # wavelength currently fixed at 510nm
        psf_vec = Cov.matvec(point_list)
        psf_vol = np.reshape(psf_vec, Cov.vol_shape)

    else:
        raise ValueError("'model_type' must be geometric or wave")

    return psf_vol

def get_Sparrow_vols(vol_points_list,calibration_list,output_file,demagnified_pitch_size=6.25):
    """
    Given a list of pairs of continuous (y,x,z) locations for two points (p_1,p_2) in the volume 
    call these 'vol_points' and the list 'vol_points_list', and a list of calibration files corresponding
    to different supersampling factors, construct the covariance matrix for each rayspread/wavespread database 
    in the list, apply it to each of the points to obtain their PSFs in the volume and save the results to HDF5.
    """
    # main HDF5 file
    try:
        print "Creating new HDF5 file:", output_file
        volumes = h5py.File(output_file,'w-')
    except:
        print "Opening existing HDF5 file:", output_file
        volumes = h5py.File(output_file,'r+')

    # main loop over points in volume and supersampling factors
    for i in xrange(len(vol_points_list)):
        
        # get volume points and create HDF5 group to store results for this loop
        vol_points = vol_points_list[i]
        vols_by_points = volumes.create_group('vol_points_'+str(i))

        # loop over calibration files for different supersampling factors
        for calibration_file in calibration_list:
            print "Analyzing:", calibration_file
            
            # create subgroup to write data to for this calibration file
            vols_by_sampling = vols_by_points.create_group('supersampling_factor_'+calibration_file.split('/')[-1].split('.')[0].split('_')[1])

            # get covariance operator
            Cov, raydb = get_Cov_from_calibration( calibration_file )

            # get point coordinates in discretized volume
            vol_coords = []
            vol_coords.append( get_voxel_coords(vol_points[0], raydb, pitch=demagnified_pitch_size))
            vol_coords.append( get_voxel_coords(vol_points[1], raydb, pitch=demagnified_pitch_size))
            print "Volume points:", vol_points
            print "Volume coordinates:",vol_coords

            # generate two psfs for vol_points and add them to get a volume containing both
            psf0 = get_psf_vol(vol_coords[0],Cov,raydb=raydb)
            psf1 = get_psf_vol(vol_coords[1],Cov,raydb=raydb)
            vol_vec = psf0 + psf1
            vol = np.reshape(vol_vec, Cov.vol_shape)
            dset = vols_by_sampling.create_dataset('Sparrow_volume', data=vol)
    volumes.close()
    return True

def get_voxel_coords(vol_points, raydb, pitch):
    """
    Get discrete voxel coordinates from continuous point sources (in microns). 
    """
    supersample_factor = raydb.supersample_factor
    x_coords = np.linspace(0,raydb.nx,num=raydb.nx+1)
    y_coords = np.linspace(0,raydb.ny,num=raydb.ny+1)
    z_coords = np.array(raydb.z_coords) 

    y = np.where(vol_points[0]<y_coords*(pitch/supersample_factor))[0][0]
    x = np.where(vol_points[1]<x_coords*(pitch/supersample_factor))[0][0]
    z = np.where(vol_points[2]<np.asarray(z_coords))[0][0]
    return (y,x,z)

def plot_Sparrow_vols(vols, vol_slice = "xy"):
        # get slice from volume
        if vol_slice == "xy":
            out_slice = vol[:,:,psf[0][2]]
        elif vol_slice == "xz":
            out_slice = vol[psf[0][0],:,:]
        elif vol_slice == "yz":
            out_slice = vol[:,psf[0][1],:]

def get_Cov_from_calibration( calibration_file ):
    """
    Given a calibration file, generate a linear covariance operator Cov = A'A.
    """
    lfcal = LightFieldCalibration.load(calibration_file)
    raydb = lfcal.rayspread_db
    lfproj = LightFieldProjection(raydb, disable_gpu = False, gpu_id = 4)
    A_op = LightFieldOperator(lfproj, raydb)
    CLO = CovarianceLinearOperator(A_op, vol_shape=(raydb.ny, raydb.nx, raydb.nz))
    Cov = LinearOperator( (raydb.nvoxels, raydb.nvoxels), matvec=CLO.matvec, rmatvec=CLO.rmatvec, dtype=np.float32 )
    Cov.vol_shape = CLO.vol_shape
    return Cov, raydb

def FFT_3D_CUDA( vol ):
    """
    Get the 3D FFT of a volume using scipy.cuda.
    """
    nx = vol.shape[1]
    ny = vol.shape[0]
    nz = vol.shape[2]
    vol.astype(np.float32)
    vol_gpu = gpuarray.to_gpu(vol) 
    F_vol_gpu = gpuarray.empty((ny, nx/2+1, nz), np.complex64) 
    plan_forward = cu_fft.Plan(vol_gpu.shape, np.float32, np.complex64) 
    cu_fft.fft(vol_gpu, F_vol_gpu, plan_forward)
    F_vol = F_vol_gpu.get()
    print 'Success status:', np.allclose(x, x_gpu.get(), atol=1e-6)
    return F_vol

def IFFT_3D_CUDA( vol_gpu, F_vol_gpu ):
    """
    Get the 3D inverse FFT of a volume using scipy.cuda.
    """
    vol_gpu_out = gpuarray.empty_like(vol_gpu) 
    plan_inverse = cu_fft.Plan(vol_gpu_out.shape, np.complex64, np.float32) 
    cu_fft.ifft(F_vol_gpu, vol_gpu_out, plan_inverse, True)
    vol_out = vol_gpu_out.get()
    print 'Success status:', np.allclose(vol_out, vol_gpu_out.get(), atol=1e-6)
    return vol_out

def test_Sparrow():
    from lflib.volume import LightFieldProjection
    from lflib.calibration import LightFieldCalibration

    print "Loading calibration data..."
    calibration_file = '/lfdata/Results/LFzfish/20121018_aTubHS_dob20121004/1.2_cyanLED15_4Hz/calibration.lfc'
    lfcal = LightFieldCalibration.load(calibration_file)
#    raydb = lfcal.rayspread_db

    # create LensletArray -- HACK, should be done through lfcal
    lenslet_array = LensletArray(lfcal.nu,lfcal.nv,lfcal.ns,lfcal.nt,lfcal.array_pitch,lfcal.magnification,lfcal.na,lfcal.medium_index,lfcal.sample_index,lfcal.lenslet_fill_factor)
    lenslet_array.objective_magnification = lenslet_array.magnification
    lenslet_array.objective_na = lenslet_array.na
    lenslet_array.ulens_pitch = lenslet_array.pitch
    lenslet_array.ulens_focal_length = lfcal.focal_length
    lenslet_array.ulens_focal_distance = lfcal.focal_distance

    print "Getting light field projection..."
    lfproj = LightFieldProjection(raydb, disable_gpu = False, gpu_id = 4)

    print "Constructing A operator..."
    A_op = LightFieldOperator(lfproj, raydb)

    print "Constructing Cov operator..."
#    CLO = CovarianceLinearOperator(A_op, vol_shape=(raydb.ny, raydb.nx, raydb.nz))
    CLO = CovarianceLinearOperator(A_op, cov_type="wave",vol_shape=(raydb.ny, raydb.nx, raydb.nz))
    Cov = LinearOperator( (raydb.nvoxels, 7), matvec=CLO.matvec, rmatvec=CLO.rmatvec, dtype=np.float32 )
    Cov.vol_shape = CLO.vol_shape

    print "Getting a PSF..."
    coords = (np.floor(raydb.ny/2),np.floor(raydb.nx/2),np.floor(raydb.nz/2))
    psf = get_psf_vol(coords, Cov, raydb=raydb, lenslet_array=lenslet_array, model_type="wave")
#    psf = get_psf_vol(coords, Cov, raydb=raydb)
    np.savez('/lfdata/Results/operator_analysis/psf.npz',psf=psf)    
    1/0
    print "Getting Sparrow criterion volumes..."
    calibration_list = os.listdir('/lfdata/Results/operator_analysis/calibration_files')
    for i in xrange(len(calibration_list)):
        calibration_list[i] = '/lfdata/Results/operator_analysis/calibration_files/'+calibration_list[i]

    # create list of points in volume for Sparrow criterion
    vol_points_list = []
    vol_points_base = [[250.0,250.0,40.0], [250.0,250.0,40.0]]
    for i in xrange(20):
        vol_points_base[1][0] +=1
        vol_points_list.append(vol_points_base) 

    # for testing
#    vol_points_list = [ [[248.0,250.0,40.0], [252.0,250.0,40.0]]] # for testing 
#    vol_points_list = [ [[248.0,250.0,40.0], [252.0,250.0,40.0]], [[247.0,250.0,40.0], [253.0,250.0,40.0]],[[246.0,250.0,40.0], [254.0,250.0,40.0]], [ [245.0,250.0,40.0], [255.0,250.0,40.0]] ]

    vols = get_Sparrow_vols(vol_points_list,calibration_list,'/lfdata/Results/operator_analysis/big_test.h5')
    1/0

#-------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    import pylab as pl
    p1s = 10**np.linspace(-10,-3,num=100)
    logLs = []; logLds = []
    for p1 in p1s:
        ks,pval,logL,logL_discrete = noisy_sensor_resolution((0,0,5e-4),(0,p1,5e-4),error_rate=1)
        logLs.append(logL)
        logLds.append(logL_discrete)
    pl.clf()
    pl.semilogx(p1s*1000, np.array(logLs)/np.max(np.array(logLs)),'r-')
    pl.semilogx(p1s*1000, np.array(logLds)/np.max(np.array(logLds)),'b--')
    pl.axhline(0,color='black')
    pl.xlabel('Distance (microns)'); #pl.xscale('log')
    pl.ylabel('Log Likelihood Ratio'); #pl.yscale('log')
    pl.grid(True)
    pl.savefig("/home/logan/Documents/Results/Resolution/sensor/logLik.png")
