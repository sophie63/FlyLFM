clear all

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','MultiSelect','on');
files=strcat(PathName,FileName);

for j=1:size(files,2)
clear D Data S Data2 out
D=MRIread(files{j});
file=files{j};
Data=D.vol;
S=size(Data);

for i=1:S(1)
    for k=1:S(3)
        parfor j=1:S(2)
            Data2(i,j,k,:)=Data(i,j,k,:)-min(Data(i,j,k,:));
        end
    end
    i
end

%Data2(Data2==NaN)=0;
out.vol=Data2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'min.nii'));


end