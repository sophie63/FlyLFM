% This script transforms a list of components stack into a colored stack

clear all

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);


parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
S=std(R(i,:));
R2=R(i,:);
R2(abs(R2)<(3*S))=0;
D2(:,:,:,i)=reshape(R2,[S1(1),S1(2),S1(3)]);
end
out.vol = D2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'thresh3std.nii'));

Dcolor=zeros(S1(1),S1(2),S1(3),3);
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

D1m=Montage3(Dcolor(:,:,:,1));
D2m=Montage3(Dcolor(:,:,:,2));
D3m=Montage3(Dcolor(:,:,:,3));
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); 
M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'coloredSumMap.PNG'));
imwrite(Dm/M, fullFileName);



out.vol = Dcolor;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'color.nii'));
