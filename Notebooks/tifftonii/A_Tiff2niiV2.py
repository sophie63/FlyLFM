from skimage import io
import nibabel as nib
import os
import numpy as np
import scipy.io as sio
import scipy.optimize
import matplotlib.pyplot as plt
import subprocess

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


## Import data from cluster
path = '/media/sophie/New Volume/100664ss1'
Dataname = path.split('/')[-1]

# Open one image to get the shape
tt = io.imread(path+'/'+Dataname+'-0002.tif') 
S = tt.shape

print ('Running...')
file_name_list = os.listdir(path)
data = np.zeros([S[0],S[1],S[2],len(file_name_list)])

for j in range(len(file_name_list)):       
    #tt = io.imread(path+'/'+Dataname+'-'+str(i+1).zfill(5)+'.tif') 
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
		plt.savefig('/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'Fit_light_on.png')
		#print ('Light on fitting plot saved at ', '/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'Fit_light_on.png')

	elif liston[-1] != len(M)-1: # use len(M)-1 if zero indexing
		print ('Estimating the time when light is off...')
		Ms=M[range(liston[len(liston)-1]-6,liston[len(liston)-1]+6)]
		res = scipy.optimize.minimize(Sq,x0=[6,3,8,-1])
		OFF=liston[len(liston)-1]-6+res.x[0]
		OFFint=np.int(np.floor(OFF))
		print("Light is off at integer time: ",OFFint)
		plt.plot(np.squeeze(Ms),'+')
		plt.plot(np.arange(0,len(Ms),0.1),[model(i,res.x[0],res.x[1],res.x[2],res.x[3]) for i in np.arange(0,len(Ms),0.1)])
		plt.title('Fit_light_off')
		plt.savefig('/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'Fit_light_off.png')
		print ('Light off fitting plot saved at ', '/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'Fit_light_off.png')

else:
	print('No light on detected.')

# Open image times
print('Opening time file...')
for f in os.listdir('/media/sophie/New Volume/Metadata/'):
    if Dataname[:6] in f:
        if f.endswith('csv'):
            TimeFile='/media/sophie/New Volume/Metadata/'+Dataname[:6]+'_.csv'
            Listfile = open(TimeFile, 'r')
            ListTime = [line.split('\n')[0] for line in Listfile.readlines()]
            Timespl=[float(ListTime[i].split(',')[2]) for i in range(1,len(ListTime))]
        elif "Original" or "Info" in f:
            TimeFile='/media/sophie/New Volume/Metadata/'+f
            with open(TimeFile, 'r') as metafile:
                lines = metafile.readlines()
            time_from_start_list = []
            for line in lines:
                if "Time_From_Start" in line:
                    split_line = line.replace('\t', ' ').replace(' = ', ' ').strip().split(' ')
                    time_from_start_list.append((int(split_line[1]),float(split_line[3])))
            Timespl = list(zip(*sorted(time_from_start_list))[1])



# Get times corresponding to images during light on (excitation light completely on : t=0)
print len(Timespl)
print ONint, OFFint
TimeOn=[Timespl[i] for i in range(ONint,(OFFint+1))]
Tinit=(ON-(ONint-1))*(Timespl[ONint]-Timespl[ONint-1])+Timespl[ONint-1]
if OFFint == len(M)-1:
    Toff = Timespl[OFFint]
else:
    Toff=(OFFint+1-OFF)*(Timespl[OFFint+1]-Timespl[OFFint])+Timespl[OFFint]

print ('Saving Time on files...')
TimeOnFinal=np.array(TimeOn)-Tinit
sio.savemat('/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'TimeFluoOn.mat', {'TimeFluoOn':TimeOnFinal})

TotalTimeOn=Toff-Tinit
sio.savemat('/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'TotalTimeOn.mat', {'TotalTimeOn':TotalTimeOn})

print ('Saving frames for which the excitation is on...')
# Keep only the frames for which the excitation is on and save
D4=np.transpose(data[:,:,:,range(ONint,(OFFint+1))],(2,1,0,3))

nim=nib.Nifti1Image(D4,np.eye(4))
fn = '/media/sophie/470fddca-e336-42c7-9d91-b91361d994ea/'+Dataname+'on.nii'
nib.save(nim, fn)
print('Finished saving file.')
print("Calling matlab function...")
os.chdir('/home/sophie/')

#subprocess.call(["matlab -nosplash -nodisplay -r \"RoughAllSteps(\'%s\')\""%(fn)], shell=True)
