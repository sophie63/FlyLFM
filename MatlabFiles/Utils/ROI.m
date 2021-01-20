clear all
%open mask
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/media');
file=strcat(PathName,FileName)
M=MRIread(file);
Mask=M.vol;

%Mask=1-Mask;%use that line if mask from segmentation editor
%Md=double(M);

%open data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/media');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;

S=size(Data);
M2=Mask./max((max(max(max(Mask)))));
SM2=size(M2);
Data(isnan(Data))=0;

M2r=reshape(M2,SM2(1)*SM2(2)*SM2(3),SM2(4));
Datar=reshape(Data,SM2(1)*SM2(2)*SM2(3),S(4));
M2rnorm=normc(M2r);

TS=M2rnorm'*Datar;

%for i=1:S(4)
    %parfor j=1:SM2(4)
        %Sm(i,j)=sum(sum(sum(M2(:,:,:,j).*squeeze(Data(:,:,:,i)))))/sum(sum(sum(M2(:,:,:,j))));
    %end
%end

save(strcat(file(1:size(file,2)-4),'TSROI.mat'),'TS')

