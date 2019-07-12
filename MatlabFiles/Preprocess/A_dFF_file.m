
%% This script does all the preprocessing steps from the 4D light reconstructed data set, to a dataset ready to be analysed using PCA and ICA
% The steps are: 
% * Detrending the data with the signal boxed averaged over 10 s. This allows to center the data and to remove decrease in fluorescence from fluorophore bleaching
% * Projecting the images along the psf depth to decrease dimentionality a
%  and concentrate the information
% * Denoising using a Kalman Filter 
% Note that files for each steps are saved so make sure that enough space
% is available on disk


clear

%% Fill up parameters here

prompt = 'What is the frame rate?';
Fr = input(prompt)

prompt = 'What is the sign of relation from deltaF/F to underlying change?';
Sdff = input(prompt)

%choose the reconstructed 4D nifti file
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);
clear D
%first detrend over 30 sec  
Unbleached_data = Sdff*dFF(Data,20*Fr);
clear Data

out.vol=Unbleached_data(:,:,:,2:(S(4)-1));
file2=strcat(file(1:size(file,2)-4),'dFF',num2str(20*Fr),'points.nii');
err = MRIwrite(out,file2);
