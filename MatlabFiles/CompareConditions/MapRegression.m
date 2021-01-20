[FileName,PathName] = uigetfile('*.nii','Select the registered template file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Masks=D.vol;
Sm=size(Masks);


Masks2=zeros(Sm(1),Sm(2),Sm(3),87);


for j=1:87
Masks2(:,:,:,j)=(Masks==j);
end
load('87to75.mat')

R2stack=zeros(Sm(1),Sm(2),Sm(3));
Coeffstack=zeros(Sm(1),Sm(2),Sm(3));

for j=1:75
R2stack=R2stack+Masks2(:,:,:,VarName3(j))*R2(j);
%Coeffstack=Coeffstack+Masks2(:,:,:,VarName3(j))*Coeff(j);
end

R2stackby2=BoxavFunction3D(R2stack);
Coeffstackby2=BoxavFunction3D(Coeffstack);
out.vol=R2stackby2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'GadR2Walkforced.nii'));
%out.vol=Coeffstackby2;
%err = MRIwrite(out,strcat(file(1:size(file,2)-4),'VglutCoeffWalk.nii'));