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
R2(abs(R2)<(3*S))=0;
D2(:,:,:,i)=reshape(R2,[S1(1),S1(2),S1(3)]);
end

 C(:,1)=[1,0,0];
 C(:,2)=[0,1,0];
 C(:,3)=[0,0,1];
 C(:,4)=[1,1,0];
 C(:,5)=[0,1,1];
 C(:,6)=[1,0,1];
 
for i=1:S1(4)
    Dcolor(:,:,:,1)=Dcolor(:,:,:,1)+D2(:,:,:,i)*C(1,(i-ceil((i-6)/6)*6));
    Dcolor(:,:,:,2)=Dcolor(:,:,:,2)+D2(:,:,:,i)*C(2,(i-ceil((i-6)/6)*6));
    Dcolor(:,:,:,3)=Dcolor(:,:,:,3)+D2(:,:,:,i)*C(3,(i-ceil((i-6)/6)*6));
end

out.vol = Dcolor;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'color.nii'));
