clear all
%open mask
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
M=MRIread(file);
Mask=M.vol;

%Mask=1-Mask;%use that line if mask from segmentation editor
%Md=double(M);

%open data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;

S=size(Data);
M2=Mask./(max(max(max(Mask))));

for i=1:S(4)
    TS(i)=sum(sum(sum(M2.*squeeze(Data(:,:,:,i)))));
end

plot(TS)
