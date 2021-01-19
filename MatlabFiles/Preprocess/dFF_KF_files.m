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

% prompt = 'What is the frame rate?';
% Fr = input(prompt)
% 
% prompt = 'What is the sign of relation from deltaF/F to underlying change?';
% Sdff = input(prompt)
Fr=200;
Sdff=1;
% prompt = 'What is the position of the focal plane in the stack?';
% z1 = input(prompt)
% 
% prompt = 'What is the distance betweeen stacks?';
% dz = input(prompt)

%% Open mask
% [FileName,PathName] = uigetfile('*.nii','Select the mask file','/home/sophie/Desktop/');
% file2=strcat(PathName,FileName)
% M=MRIread(file2);
% Mask=M.vol;
% M2=Mask./(max(max(max(Mask))));

%choose the reconstructed 4D nifti file
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','MultiSelect','on');
files=strcat(PathName,FileName);

for j=1:size(files,2)

D=MRIread(files{j});
file=files{j};
Data=D.vol;
S=size(Data);
clear D
Dtemp=mean(Data,4);
out2.vol=Dtemp;
err = MRIwrite(out2,strcat(file(1:size(file,2)-4),'temp.nii'));
%first detrend over 20 sec  
Unbleached_data = Sdff*dFF(Data,20*Fr);
clear Data

out.vol=Unbleached_data(:,:,:,2:(S(4)-1));
file2=strcat(file(1:size(file,2)-4),'dFF4000points.nii');
err = MRIwrite(out,file2);

clear out
parfor i=1:S(3)
C=squeeze(Unbleached_data(:,:,i,:));
k50=Kalman_Stack_Filter(C);
Dkf(:,:,i,:)=k50;
i
end

out.vol=Dkf(:,:,:,2:S(4)-1);
err = MRIwrite(out,strcat(file2(1:size(file2,2)-4),'kf.nii'));
clear Dkf

%S=size(Unbleached_data);
%DM=Unbleached_data;
%Dtemp=DM(:,:,:,1);

% for i=1:S(4)
% Z0=1:size(DM,3);
% Zinit=((Z0-z1)*dz+dz/2);
% z_psf=abs(Zinit)*0.239+5.46;

% % Average the stack layers when the sampling is below psf half width
% j=1;
% k=1;
% while k<=(size(DM,3))
%     if (z_psf(k)>dz)
%         nz=int8(z_psf(k)/(dz));
%         Dpsf2(:,:,j,i)=mean(DM(:,:,k:min((k+nz),size(DM,3)),i),3);
%         if i==1
%             Dtemppsf(:,:,j)=mean(Dtemp(:,:,k:min((k+nz),size(DM,3))),3);
%         end
%         Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
%         j=j+1;
%         k=k+nz+1;
%     else
%         Dpsf2(:,:,j,i)=DM(:,:,k,i);
%         if i==1
%             Dtemppsf(:,:,j)=Dtemp(:,:,k);
%         end
%         Znew(j)=Zinit(k);
%         j=j+1;
%         k=k+1;
%     end
% end
% 
% end
% clear DM
% 
% out.vol=Dpsf2;
% err = MRIwrite(out,strcat(file2(1:size(file2,2)-4),'psf.nii'));
% clear out
% 
% out2.vol=Dtemppsf;
% err = MRIwrite(out2,strcat(file2(1:size(file2,2)-4),'temp.nii'));
% clear out2

% S=size(Dpsf2);
% 
% % Use a Kalman filter to denoise the data
% parfor i=1:S(3)
% C=squeeze(Dpsf2(:,:,i,:));
% k50=Kalman_Stack_Filter(C);
% Dkf(:,:,i,:)=k50;
% i
% end
% clear Dpsf2
% 
% out.vol=Dkf(:,:,:,2:S(4)-1);
% err = MRIwrite(out,strcat(file2(1:size(file2,2)-4),'psfkf.nii'));
clear Dkf out S
clear D Data S Dkf k25 C out Dtemp Dtemppsf
end
