from skimage import io
import nibabel as nib
import os
import numpy as np
import scipy.io as sio
import scipy.optimize
#import matplotlib.pyplot as plt
import subprocess
import sys
import cv2

saving_path = "/media/NAS/Sophie/RawData/toPreprocess"
path_list = os.listdir("/media/NAS/Sophie/RawData/ss1/")

for path in path_list:
    Dataname=path.split('/')[-1]
    print ('Running on '+Dataname +" "+ "/media/NAS/Sophie/RawData/ss1/"+path+'/'+Dataname+'-00010.tif')

    # Open one image to get the shape
    tt = cv2.imread("/media/NAS/Sophie/RawData/ss1/"+path+'/'+Dataname+'-00010.tif') 
    S=tt.shape
    print(S)
    file_name_list = os.listdir("/media/NAS/Sophie/RawData/ss1/"+path)
    data = np.zeros([S[0],S[1],S[2],len(file_name_list)])
    k=0
    for j in range(len(file_name_list)):
        if os.path.exists(path+'/'+Dataname+'-'+str(j).zfill(5)+'.tif'):
            tt = io.imread(path+'/'+Dataname+'-'+str(j).zfill(5)+'.tif')
            k=k+1
            data[:,:,:,k] = tt[:][:][:]
    print ('Saving frames')
    D4=np.transpose(data,(2,1,0,3))
    nim=nib.Nifti1Image(D4,np.eye(4))
    fn = saving_path+Dataname[:6]+'/'+Dataname+'.nii'
    nib.save(nim,fn)
    print ('Finished saving individual files.')

