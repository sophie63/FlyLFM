clear all

[FileName,PathName] = uigetfile('*.nii','Select the registered template file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Masks=D.vol;
Sm=size(Masks);

[FileName,PathName] = uigetfile('*.nii','Select the average of the data Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
data=D.vol;
S=size(data);

clear D;



times = sum(sum(sum(data.*Masks)))./sum(sum(sum(Masks)));











