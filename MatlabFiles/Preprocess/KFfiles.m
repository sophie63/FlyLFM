clear

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','MultiSelect','on');
files=strcat(PathName,FileName)

for j=1:size(files,2)
clear D Data S Dkf k25 C out
D=MRIread(files{j});

Data2=D.vol;
S2=size(Data2);
Data=Data2(:,:,:,2:S2(4)-1);

S=size(Data);

parfor i=1:S(3)
C=squeeze(Data(:,:,i,:));
k75=Kalman_Stack_Filter(C);
Dkf(:,:,i,:)=k75;
i
end

out.vol=Dkf(:,:,:,2:S(4)-1);
err = MRIwrite(out,strcat(files{j}(1:size(files{j},2)-4),'kf.nii'));
end