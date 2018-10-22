clear all

% Prompt window to select target data file. Read in data file.
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol); 
clear B

% Reshape as a space*time 2D matrix 
S1=size(D);
parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end

% Demean
S2=size(R);
parfor i=1:S2(2)
    Rm(i)=mean(R(:,i));
    R(:,i)=R(:,i)-Rm(i);
end


[u,s,v]=nets_svds(R,200);
plot(cumsum(diag(s).*diag(s))/sum(diag(s).*diag(s)))

% Save PCA maps
DE=v';
for i=1:200
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'rawPCAMaps.nii'));

TS=u*s;
save(strcat(file(1:size(file,2)-4),'rawPCATS'),'TS')

