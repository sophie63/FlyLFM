clear all

%open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=B.vol;
%Masking of non brain data
%get the data before unbleaching

%reshape

S1=size(D);

parfor i=1:S1(4)
R(:,i)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end

P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=double(R).*(10000/P);

% demean
S2=size(R1);
parfor i=1:S2(1)
    R2(i,:)=R1(i,:)-mean(R1(i,:));
    %R2(i,:)=medfilt3(R2(i,:),5);
end

%need timeXspace
R=R2';

% %%Smith lines
% [uu,ss,vv]=nets_svds(R,30); % initial SVD to the top 30 components (arbitrary number fixed in MELODIC)
% vv(abs(vv)<2.3*std(vv(:)))=0;  % zero small values in PCA maps  %what is that doing here?
% stddevs=max(std(R-uu*ss*vv'),0.1);  % subtract main parts of top components from data to get normalisation
% R=R./repmat(stddevs,size(R,1),1);  % var-norm
% 
% 'hello'

%SVD
Npc=2000;
[u,s,v]=nets_svds(R,Npc);

'hello2'

%pseudo time series for fsl
DE=s*v';
for i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end

% save ICA maps
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file,'sv2000novarnorm.nii'));

save(strcat(file,'u'),'u');