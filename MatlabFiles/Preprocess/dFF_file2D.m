
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

S=size(Data)
t=1:S(3);

for (i=1:S(1))
    i
           %parfor (j=1:S(2))
           for (j=1:S(2))
        D=squeeze(Data(i,j,:));
        C=smooth(D,20*Fr);
        Cm=mean(C);
        Unbleached_data(i,j,:)=(D-C)/max([Cm,1]);
        end

    end


clear Data

out.vol=reshape(Unbleached_data(:,:,2:(S(3)-1)),[S(1) S(2) 1 S(3)-2]);
file2=strcat(file(1:size(file,2)-4),'dFF',num2str(20*Fr),'points.nii');
err = MRIwrite(out,file2);
k50=Kalman_Stack_Filter(Unbleached_data(:,:,2:(S(3)-1)));
out.vol=reshape(k50,[S(1) S(2) 1 S(3)-2]);
file2=strcat(file(1:size(file,2)-4),'dFF',num2str(20*Fr),'pointsKF.nii');
err = MRIwrite(out,file2);

