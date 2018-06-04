% This script transforms a list of components stack into a colored stack

clear all

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);

Dcolor=zeros(S1(1),S1(2),S1(3),3);
parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
S=std(R(i,:));
R2=R(i,:);
R2(abs(R2)<(2*S))=0;
D2(:,:,:,i)=reshape(R2,[S1(1),S1(2),S1(3)]);
end

out.vol = D2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'thresh2sig.nii'));