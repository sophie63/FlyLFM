from skimage import io
import nibabel as nib
import os
import numpy as np
import scipy.io as sio
import scipy.optimize
import matplotlib.pyplot as plt
import subprocess
import sys

saving_path = sys.argv[2]

path_list = [sys.argv[1]]

Fr = 100
Sdff = 1
z1 = 15
dz = 6 

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
series_name = path_list[0].split('/')[-1][:6]+'series'


# create the series dir
if not os.path.exists(saving_path):
    os.makedirs(saving_path)


# Import data from cluster
for path in path_list:
    Dataname=path.split('/')[-1]

    #inside the series dir, create folder for each dataset
    if not os.path.exists(saving_path+Dataname[:6]):
        os.makedirs(saving_path+Dataname[:6])

    # Open one image to get the shape
    tt = io.imread(path+'/'+Dataname+'-00002.tif') 
    S=tt.shape

    print ('Running on '+Dataname)
    file_name_list = os.listdir(path)
    data = np.zeros([S[0],S[1],S[2],len(file_name_list)])
    k=0

    for j in range(len(file_name_list)):
        if os.path.exists(path+'/'+Dataname+'-'+str(j).zfill(5)+'.tif'):
            tt = io.imread(path+'/'+Dataname+'-'+str(j).zfill(5)+'.tif')
            k=k+1
            data[:,:,:,k] = tt[:][:][:]
            
		
    print ('Loaded all image files.')
    print ('Estimating the time when light is on...')
    # Calculate average time series
    M=np.mean(np.mean(np.mean(data,0),0),0)
    Mav=M.mean()

    # Get approximate on and off times
    liston=[i for i in range(len(M)) if M[i]>Mav*0.8]
    liston, listoff = [],[]
    for i in range(len(M)):
        if M[i] > Mav*0.8:
            liston.append(i)
        else:
            listoff.append(i)

    ONint = 0
    ON = 0
    OFFint = len(liston)-1
    OFF = len(liston)-1

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

    		plt.savefig(saving_path+Dataname[:6]+'/'+Dataname+'_Fit_light_on.png')

    	if liston[-1] < len(M) -6 : # use len(M)-1 if zero indexing
    		print ('Estimating the time when light is off...')
    		Ms=M[range(liston[len(liston)-1]-6,liston[len(liston)-1]+6)]
    		res = scipy.optimize.minimize(Sq,x0=[6,3,8,-1])
    		OFF=liston[len(liston)-1]-6+res.x[0]
    		OFFint=np.int(np.floor(OFF))
    		print("Light is off at integer time: ",OFFint)
    		plt.plot(np.squeeze(Ms),'+')
    		plt.plot(np.arange(0,len(Ms),0.1),[model(i,res.x[0],res.x[1],res.x[2],res.x[3]) for i in np.arange(0,len(Ms),0.1)])
    		plt.title('Fit_light_off')
    		plt.savefig(saving_path+Dataname[:6]+'/'+Dataname+'_Fit_light_off.png')

    else:
    	print('No light on detected.')


#    for f in os.listdir('/media/test/Metadata/'):
#        if Dataname[:6] in f:
#            if f.endswith('csv'):
#                TimeFile='/media/test/Metadata/'+Dataname[:6]+'_.csv'
#                Listfile = open(TimeFile, 'r')
#                ListTime = [line.split('\n')[0] for line in Listfile.readlines()]
#                Timespl=[float(ListTime[i].split(',')[2]) for i in range(1,len(ListTime))]
#            elif "Original" or "Info" in f:
#                TimeFile='/media/test/Metadata/'+f
#                with open(TimeFile, 'r') as metafile:
#                    lines = metafile.readlines()
#                time_from_start_list = []
#                for line in lines:
#                    if "Time_From_Start" in line:
#                        split_line = line.replace('\t', ' ').replace(' = ', ' ').strip().split(' ')
#                        time_from_start_list.append((int(split_line[1]),float(split_line[3])))
#                #Timespl = list(zip(*sorted(time_from_start_list))[1])


    # Get times corresponding to images during light on (excitation light completely on : t=0)
    #TimeOn=[Timespl[i] for i in range(ONint,(OFFint+1))]
    #Tinit=(ON-(ONint-1))*(Timespl[ONint]-Timespl[ONint-1])+Timespl[ONint-1]
    #if OFFint == len(M)-1:
    #    Toff = Timespl[OFFint]
    #else:
     #   Toff=(OFFint+1-OFF)*(Timespl[OFFint+1]-Timespl[OFFint])+Timespl[OFFint]

    #print ('Saving Time on files...')
    #TimeOnFinal=np.array(TimeOn)-Tinit
    #sio.savemat(saving_path+Dataname[:6]+'/'+Dataname+'_TimeFluoOn.mat', {'TimeFluoOn':TimeOnFinal})
    

    #TotalTimeOn=Toff-Tinit
    #sio.savemat(saving_path+Dataname[:6]+'/'+Dataname+'_TotalTimeOn.mat', {'TotalTimeOn':TotalTimeOn})
    
    print ('Saving frames for which the excitation is on...')
    data_list.append(data[:,:,:,range(ONint,(OFFint+1))])
    # index_list has the length of each indivial file
    index_list.append(len(range(ONint,(OFFint+1))))
    # Keep only the frames for which the excitation is on and save
    D4=np.transpose(data[:,:,:,range(ONint,(OFFint+1))],(2,1,0,3))
    nim=nib.Nifti1Image(D4,np.eye(4))
    fn = saving_path+Dataname[:6]+'/'+Dataname+'on.nii'
    nib.save(nim,fn)
    print ('Finished saving individual files.')

concat = np.concatenate(data_list, axis=3)
D4_concat = np.transpose(concat, (2,1,0,3))
nim_concat=nib.Nifti1Image(D4_concat,np.eye(4))
fn = saving_path+series_name+'on.nii'
nib.save(nim_concat,fn)
print ('Finished saving concatenated file.')

# afni 3dvolreg
os.chdir(saving_path)

bashCommand = "3dvolreg -prefix "+saving_path+series_name+'reg.nii '+fn
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


#open nii image
print ('reading reg image')
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
    subprocess.call(["matlab -nosplash -nodisplay -r \"pipeline_dFF_psf_KF(\'%s\',%d,%d,%d,%d)\""%(file_name,Fr,Sdff,z1,dz)], shell=True)

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





