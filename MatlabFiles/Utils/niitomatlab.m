clear all

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
data=double(B.vol);

filer=strsplit(file,'.')
newfilename=strcat(filer{1},'.mat');

save(newfilename)
