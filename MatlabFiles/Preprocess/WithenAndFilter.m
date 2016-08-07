% This script whiten and filters the data

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

% % %simple std norm (choose between this and Smith lines below)
%   stddevs=max(std(R),0.01);  
%   R=R./repmat(stddevs,size(R,1),1);  % var-norm

%% Variance normalisation ala melodic
[uu,ss,vv]=nets_svds(R,400); % initial SVD to the top 30 components (arbitrary number fixed in MELODIC)
%vv(abs(vv)<0.4*std(vv(:)))=0;
%vv(abs(vv)<2.3*std(vv(:)))=0;
plot(diag(ss))
Rfiltered=uu(:,1:400)*vv(:,1:400)';

for i=1:S2(1)
Dpca(:,:,:,i)=reshape(Rfiltered(i,:)',[S1(1),S1(2),S1(3)]);
end
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'WFiltb.nii'));


