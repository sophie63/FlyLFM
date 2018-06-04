from skimage import io
import nibabel as nib
import os
import numpy as np
import scipy.io as sio
import scipy.optimize
import matplotlib.pyplot as plt
import subprocess


saving_path = '/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/UVvsBlue/100509series/'
path_list = ['/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/UVvsBlue/100509series/100509',
            '/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/UVvsBlue/100509series/100510']

Fr = 50
Sdff = 1
z1 = 20
dz = 6 

data_list = []
index_list = [6009, 1925]
series_name = '100509series'
fn = saving_path+series_name+'on.nii'
# create the series dir
if not os.path.exists(saving_path):
    os.makedirs(saving_path)

'''
# afni 3dvolreg
os.chdir(saving_path)

#bashCommand = "3dvolreg -prefix /media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/100630/100630ss10ntest.nii /media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/100630/100630ss1on.nii"
bashCommand = "3dvolreg -prefix "+saving_path+series_name+'reg.nii '+fn
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

'''
#open nii image
img = nib.load(saving_path+series_name+'reg.nii')
img_data = img.get_data()

# cut image based on index list, save individuals to nii and feed in to matlab dFF_psf_kf
idx_list = [0]+ list(np.cumsum(index_list))
for i in range(1,len(idx_list)):
    img = nib.Nifti1Image(img_data[:,:,:,idx_list[i-1]:idx_list[i]], np.eye(4))
    file_name = saving_path+path_list[i-1].split('/')[-1]+'reg.nii'
    
    nib.save(img, file_name)
    # change directory to open matlab
    os.chdir('/home/sophie/')
    subprocess.call(["matlab -nosplash -nodisplay -r \"pipeline_dFF_psf_KF(\'%s\',%d,%d,%d,%d)\""%(fn,Fr,Sdff,z1,dz)], shell=True)

# concatenate dkf and save
dkf_list = []
for i in range(len(path_list)):
    fn = saving_path+path_list[i].split('/')[-1]+'regdFF20spsfkf.nii'
    dkf = np.array(nib.load(fn).get_data())
    dkf_list.append(dkf)

## Something went wrong in concat
dkf_concat = np.concatenate(dkf_list, axis=3)
#dkf_D4_concat = np.transpose(dkf_concat, (2,1,0,3))
dkf_nim_concat=nib.Nifti1Image(dkf_concat,np.eye(4))
dkf_fn = saving_path+series_name+'regdFF20spsfkf_concat.nii'
nib.save(dkf_nim_concat,dkf_fn)

# pca/ica on concatenated dkf 
os.chdir('/home/sophie/')
subprocess.call(["matlab -nosplash -nodisplay -r \"pipeline(\'%s\')\""%(dkf_fn)], shell=True)





