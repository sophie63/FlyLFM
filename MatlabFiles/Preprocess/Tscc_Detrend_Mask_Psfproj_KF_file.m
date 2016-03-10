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

%Fr: frame rate
%Sdff: sign of relation from deltaF/F to underlying change
%z1: position of the focal plane in the stack
%dz:distance between z stacks
%file: nifti file
%file2: Maskfile

function []=Tscc_Detrend_Mask_Psfproj_KF_file(Fr,Sdff,z1,dz,file,file2)


M=MRIread(file2);
Mask=M.vol;
M2=Mask./(max(max(max(Mask))));

D=MRIread(file);
Data=D.vol;
S=size(Data);
clear D

%first detrend over 7sec 
Unbleached_data = Sdff*Detrend(Data,Fr);
clear Data

out.vol=Unbleached_data(:,:,:,2:(S(4)-1));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U7s.nii'));

clear out

%%Mask the detrended data
S=size(Unbleached_data);

parfor i=1:S(4)
    DM(:,:,:,i)=M2.*Unbleached_data(:,:,:,i);
end
clear Unbleached_data

out.vol=DM;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U7sM.nii'));
clear out

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
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j,i)=DM(:,:,k,i);
        Znew(j)=Zinit(k);
        j=j+1;
        k=k+1;
    end
end

i
end
clear DM

out.vol=Dpsf2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U7sMpsf.nii'));
clear out

S=size(Dpsf2);

% Use a Kalman filter to denoise the data
parfor i=1:S(3)
C=squeeze(Dpsf2(:,:,i,:));
k50=Kalman_Stack_Filter(C);
Dkf(:,:,i,:)=k50;
i
end
clear Dpsf2

out.vol=Dkf(:,:,:,2:S(4)-1);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'U7sMpsfkf.nii'));



