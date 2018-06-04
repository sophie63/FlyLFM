clear

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','MultiSelect','on');
files=strcat(PathName,FileName);


D=MRIread(files{1});
Data=D.vol;
S=size(Data);

for j=2:size(files,2)
clear D1 
D1=MRIread(files{j});
Data1=D1.vol;
S=size(Data);

Data=cat(4,Data,Data1);

end

% save
out.vol=Data;
err = MRIwrite(out,strcat(files{1}(1:size(files{1},2)-4),'concat.nii'));
