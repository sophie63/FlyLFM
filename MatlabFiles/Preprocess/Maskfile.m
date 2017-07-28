clear all
%open mask
[FileName,PathName] = uigetfile('*.nii','Select the Mask file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
M=MRIread(file);
Mask=M.vol;
M2=Mask./(max(max(max(Mask))));
%M2=1-M2;%use that line if mask from segmentation editor
%Md=double(M);


[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)

D=MRIread(file);
Data=D.vol;

S=size(Data);


parfor i=1:S(4)
    DM(:,:,:,i)=M2.*Data(:,:,:,i);
end

out.vol=DM;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'M.nii'));

