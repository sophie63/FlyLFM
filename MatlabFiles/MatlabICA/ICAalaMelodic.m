% This script performs PCA and spatial ICA in a manner similar to melodic from FSL 

clear all

% Open Data
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
NPCfilt=60;
[uu,ss,vv]=nets_svds(R,NPCfilt); 
% initial SVD to the top components
vv(abs(vv)<0.4*std(vv(:)))=0;
%vv(abs(vv)<2.3*std(vv(:)))=0;
%vv(abs(vv)<0.023*std(vv(:)))=0;
stddevs=max(std(R-uu*ss*vv'),1);  
% Subtract main parts of top components from data to get normalisation
% Which leaves most main components not normalized
% Note also that small values are not increased (normalized by one)
R=R./repmat(stddevs,size(R,1),1);  % Var-norm

% Check dimentionality
[uu,ss,vv]=nets_svds(R,1000); 
plot(log(diag(ss)))
% pause

%% SVD
% Use twice the inflexion point 
Npc=150;
[u,s,v]=nets_svds(R,Npc);

% Save PCA maps
DE=s*v';
for i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCAMaps.nii'));

save(strcat(file(1:size(file,2)-4),'PCATS'),'u')

%pause
% Check the PCs then, identify the components that look like movement and
% noise and remove them before ICA
Id=[1,3,5];
goodPC=[2,4,6:150];
Npc=size(goodPC,2);
newv=v(:,goodPC);

%% ICA
[icasig, A, W] = fastica (newv','approach','symm','epsilon', 0.001);

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

TS=u(:,goodPC)*s(goodPC,goodPC)*A;
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
err = MRIwrite(out,strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'IC.nii'));

save(strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'TS'),'TSo')

% Next step is opening the maps and time series in an ipython notebook for
% manual sorting
% Zscore maps
for i=1:Npc
    GM1vn(:,i)=GM(:,i)/(sqrt(var(GM(:,i))));
end

GMz=GM;
GMz(GM1vn<2.5)=0;

% Remove bad PCs
Rnew=R-u(:,Id)*s(Id,Id)*v(:,Id)';
for j=1:Npc
parfor i=1:S1(4)
TSzmap(i,j)=mean(Rnew(i,:).*GMz(:,j)');
end
j
end

TSzmapo=TSzmap(:,Order(Npc:-1:1));
save(strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'TSzmap'),'TSzmapo')

figure
hold on
for i=1:Npc
    plot(TSo(:,i)/sqrt(var(TSo(:,i)))+i*10,'r')
    plot(TSzmapo(:,i)/sqrt(var(TSzmapo(:,i)))+i*10+5,'b')
end
