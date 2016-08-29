%% This script does all the preprocessing steps from the 4D light reconstructed data set, to a dataset ready to be analysed using PCA and ICA
% The steps are: 
% * Detrending the data with the signal boxed averaged over 7 s. This allows to center the data and to remove decrease in fluorescence from fluorophore bleaching
% * Masking the data using a thresholded version of an image (prepared in
% imageJ)
% * Projecting the images along the psf depth to decrease dimentionality a
%  and concentrate the information
% * Denoising using a Kalman Filter 
% Note that files for each steps are saved so make sure that enough space
% is available on disk


clear all

prompt = 'What is the frame rate?';
Fr = input(prompt)

prompt = 'What is the sign of relation from deltaF/F to underlying change?';
Sdff = input(prompt)

prompt = 'What is the position of the focal plane in the stack?';
z1 = input(prompt)

prompt = 'What is the distance betweeen stacks?';
dz = input(prompt)


%open mask
[FileName,PathName] = uigetfile('*.nii','Select the mask file','/home/sophie/Desktop/');
file2=strcat(PathName,FileName)
M=MRIread(file2);
Mask=M.vol;
M2=Mask./(max(max(max(Mask))));

%choose the reconstructed 4D nifti file
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);
clear D


% First detrend over 10 sec 
Unbleached_data = Sdff*Detrend(Data,1/Fr);

clear Data

out.vol=Unbleached_data(:,:,:,2:(S(4)-1));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U10s.nii'));

clear out

%%Mask the detrended data
S=size(Unbleached_data);

parfor i=1:S(4)
    DM(:,:,:,i)=M2.*Unbleached_data(:,:,:,i);
end
clear Unbleached_data
% 
out.vol=DM;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U10sM.nii'));
clear out

Dtemp=DM(:,:,:,1);

for i=1:S(4)
Z0=1:size(DM,3);
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;

% Average the stack layers when the sampling is below psf half width
j=1;
k=1;
while k<=(size(DM,3))
    if (z_psf(k)>dz)
        nz=int8(z_psf(k)/(dz));
        Dpsf2(:,:,j,i)=mean(DM(:,:,k:min((k+nz),size(DM,3)),i),3);
        if i==1
            Dtemppsf(:,:,j)=mean(Dtemp(:,:,k:min((k+nz),size(DM,3))),3);
        end
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j,i)=DM(:,:,k,i);
        if i==1
            Dtemppsf(:,:,j)=Dtemp(:,:,k);
        end
        Znew(j)=Zinit(k);
        j=j+1;
        k=k+1;
    end
end
end
clear DM

out.vol=Dpsf2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U10sMpsf.nii'));
clear out

out2.vol=Dtemppsf;%   Detailed explanation goes here
err = MRIwrite(out2,strcat(file(1:size(file,2)-4),'temp.nii'));
clear out2

S=size(Dpsf2);

% Use a Kalman filter to denoise the data
parfor i=1:S(3)
C=squeeze(Dpsf2(:,:,i,:));
k50=Kalman_Stack_Filter(C);
Dkf(:,:,i,:)=k50;
end
clear Dpsf2

out.vol=Dkf(:,:,:,2:S(4)-1);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U10sMpsfkf.nii'));




