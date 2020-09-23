
% open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)

D=MRIread(file);
Data=D.vol;

S=size(Data);
d=16;
Ds=zeros(S(1),S(2),floor(S(3)/d));
for i=1:floor(S(3)/d)
Ds(:,:,i)=Data(:,:,i*d-3)+Data(:,:,i*d-2)+Data(:,:,i*d-1)+Data(:,:,i*d);
end

out.vol=Ds;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'TimeBin4.nii'));
