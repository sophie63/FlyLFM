clear

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);
clear D

Unbleached_data = - Detrend( Data,200);
clear Data

out.vol=Unbleached_data(:,:,:,2:(S(4)-1));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U7s.nii'));
