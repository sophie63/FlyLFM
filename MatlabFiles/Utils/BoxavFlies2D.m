
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','MultiSelect','on');
files=strcat(PathName,FileName);

for j=1:size(files,2)
clear D Data S Data2 out
D=MRIread(files{j});
file=files{j};
% open Data

D=MRIread(file);
Data=D.vol;

S=size(Data);

Ds=zeros(floor(S(1)/2), floor(S(2)/2),S(3));
for i=1:S(3)
        X=squeeze(Data(1:floor(S(1)/2)*2,1:floor(S(2)/2)*2,i));
        Xs=reshape(sum(sum(reshape(X, 2, floor(S(1)/2), 2, floor(S(2)/2)), 1), 3), floor(S(1)/2), floor(S(2)/2));
        Ds(:,:,i)=Xs;
    i
end


out.vol=Ds;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'boxav.nii'));
end
