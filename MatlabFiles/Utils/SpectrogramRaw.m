
clear all

% Open Data
%[FileName,PathName] = uigetfile('*.h5','Select the  file');
%file=strcat(PathName,FileName)
%B=h5read(file,'/imaging');
%D=double(B);
% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=squeeze(double(B.vol));


clear B

% Reshape as a space*time 2D matrix 
S1=size(D);
parfor i=1:S1(3)
R(i,:)=reshape(D(:,:,i),[1,S1(1)*S1(2)]);
end

M=mean(R(1000:6000,:),2);
C=smooth(M,10*100);
Cm=mean(C);
M=(M-C);

params.Fs=100;
movingwin=[5 1];
[S,t,f]=mtspecgramc(M,movingwin,params);
figure
plot_matrix(S,t,f)