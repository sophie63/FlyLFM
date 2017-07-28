% This script performs PCA and spatial ICA with a normalisation similar to melodic from FSL 

clear all

% Open Data
%[FileName,PathName] = uigetfile('*.h5','Select the  file');
%file=strcat(PathName,FileName)
%B=h5read(file,'/imaging');
%D=double(B);
% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=squeeze(double(B.vol));


clear B

% Reshape as a space*time 2D matrix 
S1=size(D);
parfor i=1:S1(3)
R(i,:)=reshape(D(:,:,i),[1,S1(1)*S1(2)]);
end
clear D
R(isnan(R))=0;
R(isinf(R))=0;
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
NPCfilt=60;
[uu,ss,vv]=nets_svds(R,NPCfilt); 
% initial SVD to the top components
vv(abs(vv)<0.4*std(vv(:)))=0;
%vv(abs(vv)<2.3*std(vv(:)))=0;
%vv(abs(vv)<0.023*std(vv(:)))=0;
stddevs=max(std(R-uu*ss*vv'),1);  
% Subtract main parts of tcop components from data to get normalisation
% Which leaves most main components not normalized
% Note also that small values are not increased (normalized by one)
R=R./repmat(stddevs,size(R,1),1);  % Var-norm

% Check dimentionality
%[uu,ss,vv]=nets_svds(R,1000); 
%plot(log(diag(ss)))
% Typically draw a line using the tail and choose Npc at the point the
% curve deviates from this (should be higher than the inflexion point)
%prompt = 'How many components do you want?';
%Npc = input(prompt)
Npc=225;

%% SVD
[u,s,v]=nets_svds(R,Npc);

% Save PCA maps
DE=s*v';
for i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2)]);
end
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCAMaps.nii'));

save(strcat(file(1:size(file,2)-4),'PCATS'),'u')

%pause
% Check the PCs then, identify the components that look like movement and
% noise and remove them before ICA
%prompt = 'What components correspond to movement? ';
%Id1 = input(prompt)
%Id2 = input(prompt)
%Id3 = input(prompt)
%Id4 = input(prompt)
%Id5 = input(prompt)
%Id6 = input(prompt)
%Id7 = input(prompt)
%Id=[Id1,Id2,Id3,Id4,Id5,Id6,Id7]

goodPC1=[1:Npc];
%goodPC1(Id)=nan;
goodPC=goodPC1(~isnan(goodPC1));
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
Dica(:,:,:,i)=reshape(GMo(:,i),[S1(1),S1(2)]);
end

% Save ICA maps and time series
out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'IC.nii'));

save(strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'TS'),'TSo')

% Use Zscored maps as region of interest
for i=1:Npc
    GM1vn(:,i)=GM(:,i)/(sqrt(var(GM(:,i))));
end

GMz=GM;
GMz(GM1vn<2.5)=0;

GM1vn(GM1vn<2.5)=0;


for j=1:Npc
parfor i=1:S1(3)
TSzmap(i,j)=mean(R0(i,:).*GMz(:,j)');
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


% Make a summary colored map of the ROIs

GMzcolor=zeros(size(GMz,1),3);
for i=1:Npc
    C=rand(3,1);
    GMzcolor(:,1)=GMzcolor(:,1)+GM1vn(:,i)*C(1);
    GMzcolor(:,2)=GMzcolor(:,2)+GM1vn(:,i)*C(2);    
    GMzcolor(:,3)=GMzcolor(:,3)+GM1vn(:,i)*C(3);    
end


% Reform volumes
for i=1:3
Dcolor(:,:,i)=reshape(GMzcolor(:,i),[S1(1),S1(2)])/max(GMzcolor(:,i));
end

% Save ICA maps and time series

fullFileName = strcat(file(1:size(file,2)-4),num2str(Npc),'Smith0_4_',num2str(NPCfilt),'ICcolored.png');
imwrite(Dcolor, fullFileName);
% Next step is opening the maps and time series in an ipython notebook for
% manual sorting
