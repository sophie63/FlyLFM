clear

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName);
D=MRIread(file);
Data=D.vol;
S=size(Data);
clear D

Npoint=10;
C=Data;
H=Data;

D=zeros(S(4));

for (i=1:S(1))
    i
    for (j=1:S(2))
        parfor(k=1:S(3))
        D=squeeze(Data(i,j,k,:));
        C(i,j,k,:)=smooth(D,Npoint);
        H(i,j,k,:)=Data(i,j,k,:)-C(i,j,k,:);
        end

    end
end


out.vol=C(:,:,:,2:(S(4)-1));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'Smoothed10.nii'));

out.vol=H(:,:,:,2:(S(4)-1));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'HighPass10.nii'));
