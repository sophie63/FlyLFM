% This script performs PCA and spatial ICA with a normalization similar to melodic from FSL 

clear all

% Prompt window to select target data file. Read in data file. Open file
% with movements already removed
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol); 
clear B

% Reshape as a space*time 2D matrix 
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

%% Simple std norm (choose between this and Smith lines below)
%   stddevs=max(std(R),0.01);  
%   R=R./repmat(stddevs,size(R,1),1);  % var-norm

%% Variance normalisation ala melodic
% NPCfilt=60;
% [uu,ss,vv]=nets_svds(R,NPCfilt); 
% % initial SVD to the top components
% vv(abs(vv)<0.4*std(vv(:)))=0;
% %vv(abs(vv)<2.3*std(vv(:)))=0;
% %vv(abs(vv)<0.023*std(vv(:)))=0;
% stddevs=max(std(R-uu*ss*vv'),1);  
% % Subtract main parts of tcop components from data to get normalisation
% % Which leaves most main components not normalized
% % Note also that small values are not increased (normalized by one)
% R=R./repmat(stddevs,size(R,1),1);  % Var-norm

% Check dimentionality
[uu,ss,vv]=nets_svds(R,1000); 
plot(log(diag(ss)))
% Choose NPC at twice the elbow (empirically found that allows to get most
% activity related components without too many noise components
Spectrum=log(diag(ss));
Diffvalue=-(max(Spectrum)-min(Spectrum))/1000;
DiffFunc=smooth(diff(smooth(Spectrum,20)));
Npc=2*find(DiffFunc>Diffvalue, 1 );


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

% pause
% Check the PCs then, identify the components that look like movement and
% noise and remove them before ICA
% prompt = 'What components correspond to movement? ';
% Id1 = input(prompt)
% Id2 = input(prompt)
% Id3 = input(prompt)
% Id4 = input(prompt)
% Id5 = input(prompt)
% Id6 = input(prompt)
% Id7 = input(prompt)
% Id=[Id1,Id2,Id3,Id4,Id5,Id6,Id7]
% 
% goodPC1=[1:Npc];
% goodPC1(Id)=nan;
% goodPC=goodPC1(~isnan(goodPC1));
% Npc=size(goodPC,2);
% newv=v(:,goodPC);

%% ICA
[icasig, A, W] = fastica (v','approach','symm','epsilon', 0.001);

% form back ica maps
GM=icasig';
GMv=GM./repmat(std(GM),size(GM,1),1);

%Correct sign so that mean of the positive side is larger than the mean of
%the negative side after zscoring

GMzp=GMv-2;
GMzp(GMzp<0)=0;
GMzpm=mean(GMzp);
GMzn=GMv+2;
GMzn(GMzn>0)=0;
GMznm=mean(GMzn);
GMs=GM.*repmat(sign(GMzpm+GMznm),size(GM,1),1); 

TS=u*s*A;
TSs=TS.*repmat(sign(GMzpm+GMznm),size(TS,1),1);

% Get back df/f
for i=1:Npc
TSmean(:,i)=TSs(:,i)/mean(GMs(:,i));
end

% Reform volumes
for i=1:Npc
Dica(:,:,:,i)=reshape(GM(:,i),[S1(1),S1(2),S1(3)]);
end

% Save ICA maps and time series
out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),num2str(Npc),'novarnorm',num2str(Npc),'IC.nii'));

save(strcat(file(1:size(file,2)-4),num2str(Npc),'novarnorm',num2str(Npc),'TS'),'TSs')

