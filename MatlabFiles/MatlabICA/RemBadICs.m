clear all

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);

load(strcat(file(1:size(file,2)-6),'TS.mat'))
TS=TSo;
IndBad=[1:12,16:31,35,37:44,47:49,51:53,56:58,60:67,70:75,77:80,82:85,87:88,90,91,93:99,102,103,105:108,111:117,119:127,129:179];
Indones=ones(1,S1(4));
Indones(IndBad)=0;
Indtot=1:S1(4);
IndGood=Indtot(logical(Indones));
TSo=TS(:,IndGood);
Dica=D(:,:,:,IndGood);

out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'goodIC.nii'));

save(strcat(file(1:size(file,2)-4),'goodTS'),'TSo')
clear TSo
TSo=TS(:,IndBad);
Dica=D(:,:,:,IndBad);
out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'badIC.nii'));

save(strcat(file(1:size(file,2)-4),'badTS'),'TSo')
