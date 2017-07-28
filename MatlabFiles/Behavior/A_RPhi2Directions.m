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

Lefta=-X;
Lefta(Lefta<0)=0;
Righta=X;
Righta(Righta<0)=0;
Straighta=-Y;

Left=zeros(1,S(4)+1);
Left(1)=Lefta(1);
Left(2:S(4))=interp1(1:S(4),squeeze(Lefta),((2:S(4))-0.5)',[]);
Left(S(4)+1)=Lefta(S(4));

Right=zeros(1,S(4)+1);
Right(1)=Righta(1);
Right(2:S(4))=interp1(1:S(4),squeeze(Righta),((2:S(4))-0.5)',[]);
Right(S(4)+1)=Righta(S(4));

Straight=zeros(1,S(4)+1);
Straight(1)=Straighta(1);
Straight(2:S(4))=interp1(1:S(4),squeeze(Straighta),((2:S(4))-0.5)',[]);
Straight(S(4)+1)=Straighta(S(4));

endfile=strcat(file(1:size(file,2)-4),'Directions.mat')
save(endfile,'Left','Right','Straight')

