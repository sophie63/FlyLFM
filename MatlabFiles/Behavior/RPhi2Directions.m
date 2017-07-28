clear all
% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B

S=size(D)
Rraw=D(:,:,1,:);
Rraw(isnan(Rraw))=0;
Phiraw=D(:,:,2,:);
Phiraw(isnan(Phiraw))=0;
Xraw=Rraw.*sin(Phiraw);
Yraw=Rraw.*cos(Phiraw);

X=squeeze(mean(mean(Xraw)));
Y=squeeze(mean(mean(Yraw)));

Left=-X;
Left(Left<0)=0;
Right=X;
Right(Right<0)=0;
Straight=-Y;
