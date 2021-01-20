from skimage import io
import nibabel as nib
import os
import numpy as np
import scipy.io as sio
import scipy.optimize
import matplotlib.pyplot as plt
import subprocess


saving_path = '/media/sophie/Samsung_T5/ToPutTogether/'

path_list = ['/media/sophie/Samsung_T5/ToPutTogether/Data244ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/Data346bss1',
             '/media/sophie/Samsung_T5/ToPutTogether/Data371ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/Data100761ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB7ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB41ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB44ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB99ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB117ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB122ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataB123ss1',
             '/media/sophie/Samsung_T5/ToPutTogether/DataTUM44ss1']


# Model for fitting onset and offset
def model(x,a,b,c,d):
    if x<a:
        return b
    elif x<c:
        return b+(x-a)*d
    else:
        return (c-a)*d+b
def Sq(X):
    return sum([(model(i,X[0],X[1],X[2],X[3])-Ms[i])**2 for i in range(len(Ms))])   


data_list = []
index_list = []

# create the series dir
if not os.path.exists(saving_path):
    os.makedirs(saving_path)


# Import data from cluster
for path in path_list:
    Dataname=path.split('/')[-1]

    # Open one image to get the shape
    tt = io.imread(path+'/'+Dataname+'-00002.tif') 
    S=tt.shape

    print ('Running on '+Dataname)
    file_name_list = os.listdir(path)
    data = np.zeros([S[0],S[1],S[2],len(file_name_list)])

    for j in range(len(file_name_list)):
        if os.path.exists(path+'/'+Dataname+'-'+str(j).zfill(5)+'.tif'):
            tt = io.imread(path+'/'+Dataname+'-'+str(j).zfill(5)+'.tif')
            data[:,:,:,j] = tt[:][:][:]

    print ('Loaded all image files.')
    print ('Estimating the time when light is on...')
    # Calculate average time series
    M=np.mean(np.mean(np.mean(data,0),0),0)
    Mav=M.mean()

    # Get approximate on and off times
    liston=[i for i in range(len(M)) if M[i]>Mav*0.7]
    liston, listoff = [],[]
    for i in range(len(M)):
        if M[i] > Mav*0.7:
            liston.append(i)
        else:
            listoff.append(i)

    ONint = 0
    ON = 0
    OFFint = len(M)-1
    OFF = len(M)-1

    if len(liston)>0:
    	if liston[0] != 0:
    		# Model onset and find precise onset time
    		Ms=M[range(liston[0]-8,liston[0]+8)]
    		res = scipy.optimize.minimize(Sq,x0=[7,0.3,9,0.7])
    		ON=liston[0]-8+res.x[2]
    		ONint=np.int(np.ceil(ON))
    		print ("Light is on at integer time: ", ONint)
    		plt.plot(np.squeeze(M[range(liston[0]-8,liston[0]+8)]),'+')
    		plt.plot(np.arange(0,len(Ms),0.1),[model(i,res.x[0],res.x[1],res.x[2],res.x[3]) for i in np.arange(0,len(Ms),0.1)])
    		plt.title('Fit_light_on')
    		plt.savefig(saving_path+Dataname+'_Fit_light_on.png')

    		
    	elif liston[-1] != len(M) -1 : # use len(M)-1 if zero indexing
    		print ('Estimating the time when light is off...')
    		Ms=M[range(liston[len(liston)-1]-6,liston[len(liston)-1]+6)]
    		res = scipy.optimize.minimize(Sq,x0=[6,3,8,-1])
    		OFF=liston[len(liston)-1]-6+res.x[0]
    		OFFint=np.int(np.floor(OFF))
    		print("Light is off at integer time: ",OFFint)
    		plt.plot(np.squeeze(Ms),'+')
    		plt.plot(np.arange(0,len(Ms),0.1),[model(i,res.x[0],res.x[1],res.x[2],res.x[3]) for i in np.arange(0,len(Ms),0.1)])
    		plt.title('Fit_light_off')
    		plt.savefig(saving_path+Dataname+'_Fit_light_off.png')

    else:
    	print('No light on detected.')

    print ('Saving frames for which the excitation is on...')

    # index_list has the length of each indivial file
    index_list.append(len(range(ONint,(OFFint+1))))
    # Keep only the frames for which the excitation is on and save
    D4=np.transpose(data[:,:,:,range(ONint,(OFFint+1))],(2,1,0,3))
    nim=nib.Nifti1Image(D4,np.eye(4))
    fn = saving_path+Dataname+'.nii'
    data_list.append(fn)
    nib.save(nim,fn)
    print ('Finished saving individual files.')







