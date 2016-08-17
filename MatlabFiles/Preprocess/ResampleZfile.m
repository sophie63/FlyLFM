<<<<<<< HEAD
clear all

% open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)

D=MRIread(file);
Data=D.vol;

S=size(Data);

for k=1:S(1)
for j=1:S(2)
parfor i=1:S(4)
    Dres(k,j,:,i)=resample(Data(k,j,:,i),2,1);
end
end
end

out.vol=Dres;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'zx2.nii'));
||||||| merged common ancestors
=======
clear all

% open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)

D=MRIread(file);
Data=D.vol;

S=size(Data);

for k=1:S(1)
for j=1:S(2)
parfor i=1:S(4)
    Dres(k,j,:,i)=resample(Data(k,j,:,i),2,1);
end
end
k
end

out.vol=Dres;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'zx2.nii'));
>>>>>>> 57759d29a73933bc92c6aaedf808d46b8fabe631
