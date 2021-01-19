%% This script does all the preprocessing steps from the 4D light reconstructed data set, to a dataset ready to be analysed using PCA and ICA
% The steps are: 
% * Detrending the data with the signal boxed averaged over 10 s. This allows to center the data and to remove decrease in fluorescence from fluorophore bleaching
% * Projecting the images along the psf depth to decrease dimentionality a
%  and concentrate the information
% * Denoising using a Kalman Filter 
% Note that files for each steps are saved so make sure that enough space
% is available on disk


clear

Fr = 200;
Sdff = 1;

[FileName,PathName] = uigetfile('*.nii','Select the Nifti files','MultiSelect','on');
files=strcat(PathName,FileName);

for j=1:size(files,2)

D=MRIread(files{j});
file=files{j};
Data=D.vol;
S=size(Data);
clear D
Template=mean(Data,4);
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


D=Dkf(3:(S(1)-2),3:(S(2)-2),3:(S(3)-2),:);

clear Dkf
% Reshape as a space*time 2D matrix 
S1=size(D);
R=reshape(D,[S1(1)*S1(2)*S1(3)],S1(4))';
clear D
R0=R;
% Demean
S2=size(R);
R=R-mean(R,1);
% Normalize
P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=R.*(10000/P);
clear R
R=R1;

%% Simple std norm (choose between this and Smith lines below)
%   stddevs=max(std(R),0.01);  
%   R=R./repmat(stddevs,size(R,1),1);  % var-norm

%% Variance normalisation ala melodic
% Check dimentionality
[u,s,v]=nets_svds(R,min(1000,S1(4))); 
plot(log(diag(s)))
% Choose NPCfilt at the elbow 
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
% Subtract main parts of top components from data to get normalisation
% Which leaves most main components not normalized
% Note also that small values are not increased (normalized by one)
R=R./repmat(stddevs,size(R,1),1);  % Var-norm

[u,s,v]=nets_svds(R,min(1000,S1(4))); 
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
%for i=1:Npc
%Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
%end
Dpca=reshape(DE',[S1(1),S1(2),S1(3),Npc]);
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

TS=u*s*A;
TSs=TS.*repmat(sign(GMzpm+GMznm),size(TS,1),1);

%reorder maps by variance 
[varo,Order]=sort(var(TS));
%plot(varo)
TSo=TSs(:,Order(Npc:-1:1));
GMo=GMs(:,Order(Npc:-1:1));


% Save ICA maps and time series
Dica=reshape(GMo,[S1(1),S1(2),S1(3),Npc]);
out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),num2str(Npc),'IC.nii'));

save(strcat(file(1:size(file,2)-4),num2str(Npc),'TS'),'TSo')
TemplICA=max(reshape(GMs,[S1(1),S1(2),S1(3),Npc]),[],4);
out.vol =TemplICA;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),num2str(Npc),'Template.nii'));

% Get ROI time series
GMv=GMo./repmat(std(GMo),size(GMo,1),1);
GMzp=GMv-2;
GMzp(GMzp<0)=0;
TSzmapo=R0*GMzp;
save(strcat(file(1:size(file,2)-4),num2str(Npc),'TSzmap'),'TSzmapo')

figure
hold on
for i=1:Npc
    plot(TSo(:,i)/sqrt(var(TSo(:,i)))+i*10,'r')
    plot(TSzmapo(:,i)/sqrt(var(TSzmapo(:,i)))+i*10+5,'b')
end

savefig(strcat(file(1:size(file,2)-4),'TimeSeries'))


clear Dkf out S TSzmapo 
clear D Data S Dkf k25 C out Dtemp Dtemppsf
end
