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

% Reshape as a space*time 2D matrix 
clear D
D=Dkf;
S1=size(D);
parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end
clear D
R0=R;

% Demean
S2=size(R);
parfor i=1:S2(2)
    Rm(i)=mean(R(:,i));
    R(:,i)=R(:,i)-Rm(i);
end

P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=R.*(10000/P);
clear R

R=R1;

%% Simple std norm (choose between this and Smith lines below)
%   stddevs=max(std(R),0.01);  
%   R=R./repmat(stddevs,size(R,1),1);  % var-norm

%% Variance normalisation ala melodic
% Check dimentionality
[u,s,v]=nets_svds(R,1000); 
plot(log(diag(s)))
% Choose NPCfilt at the elbow (empirically found that allows to get most
% activity related components without too many noise components
Spectrum=log(diag(s));
Diffvalue=-(max(Spectrum)-min(Spectrum))/1000;
DiffFunc=smooth(diff(smooth(Spectrum,20)));
NPCfilt=find(DiffFunc>Diffvalue, 1 )

uu=u(:,1:NPCfilt);
vv=v(:,1:NPCfilt);
ss=s(1:NPCfilt,1:NPCfilt);

vv(abs(vv)<0.4*std(vv(:)))=0;
%vv(abs(vv)<2.3*std(vv(:)))=0;
%vv(abs(vv)<0.023*std(vv(:)))=0;
stddevs=max(std(R-uu*ss*vv'),1);  
% Subtract main parts of tcop components from data to get normalisation
% Which leaves most main components not normalized
% Note also that small values are not increased (normalized by one)
R=R./repmat(stddevs,size(R,1),1);  % Var-norm



[u,s,v]=nets_svds(R,1000); 
plot(log(diag(s)))
% Choose NPCfilt at 2*the elbow (empirically found that allows to get most
% activity related components without too many noise components
Spectrum=log(diag(s));
Diffvalue=-(max(Spectrum)-min(Spectrum))/1000;
DiffFunc=smooth(diff(smooth(Spectrum,20)));
Npc=2*find(DiffFunc>Diffvalue, 1 )

%% SVD
[u,s,v]=nets_svds(R,Npc);

% Save PCA maps
DE=s*v';
for i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCAMaps.nii'));

save(strcat(file(1:size(file,2)-4),'PCATS'),'u')


%% ICA
[icasig, A, W] = fastica (v','approach','symm','epsilon', 0.001);


%% Correct sign so that mean of the positive side is larger than the mean of
% the negative side after zscoring
% form back ica maps

GM=icasig';
GMv=GM./repmat(std(GM),size(GM,1),1);
GMzp=GMv-2;
GMzp(GMzp<0)=0;
GMzpm=mean(GMzp);
GMzn=GMv+2;
GMzn(GMzn>0)=0;
GMznm=mean(GMzn);
GMs=GM.*repmat(sign(GMzpm+GMznm),size(GM,1),1); 

TS=u*s*A*stdev;
TSs=TS.*repmat(sign(GMzpm+GMznm),size(TS,1),1);

%reorder maps by variance 
[varo,Order]=sort(var(TS));
%plot(varo)
TSo=TSs(:,Order(Npc:-1:1));
GMo=GMs(:,Order(Npc:-1:1));

% Reform volumes
for i=1:size(goodPC,2)
Dica(:,:,:,i)=reshape(GMo(:,i),[S1(1),S1(2),S1(3)]);
end

% Save ICA maps and time series
out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),num2str(Npc),num2str(NPCfilt),'IC.nii'));

save(strcat(file(1:size(file,2)-4),num2str(Npc),num2str(NPCfilt),'TS'),'TSo')


% Next step is opening the maps and time series in an ipython notebook for
% manual sorting

clear Dkf out S
clear D Data S Dkf k25 C out Dtemp Dtemppsf
end
