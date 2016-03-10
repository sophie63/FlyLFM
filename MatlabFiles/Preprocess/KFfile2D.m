clear

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

%drop vid 

%parfor i=1:S(3)
%C=squeeze(Data(:,:,i,:));
k75=Kalman_Stack_Filter(Data);
%Dkf(:,:,i,:)=k75;
%i
%end

%out.vol=Dkf(:,:,:,2:S(4)-1);
out.vol=Dkf(:,:,:);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'kf.nii'));