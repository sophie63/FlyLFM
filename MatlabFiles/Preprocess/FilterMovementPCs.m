% This script performs PCA and spatial ICA in a manner similar to melodic from FSL 

clear all

% Open Data
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
clear D

%% SVD
% Use twice the inflexion point 
Npc=150;
[u,s,v]=nets_svds(R,Npc);

% Save PCA maps
DE=s*v';
parfor i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCAMaps.nii'));

prompt = 'What are the components that correspond to movement? ';
Id1 = input(prompt)
Id2 = input(prompt)
Id3 = input(prompt)
Id=[Id1,Id2,Id3]
%Id=1:Npc;
%Id(x)=nan;
%Id=Id(~isnan(Id));
Rnew=R-u(:,Id)*s(Id,Id)*v(:,Id)';

% Save filtered data
parfor i=1:S1(4)
Dnew(:,:,:,i)=reshape(Rnew(i,:),[S1(1),S1(2),S1(3)]);
end
out4.vol = Dnew;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCAF.nii'));


% Rnewtest=R-u*s*v';
% 
% % Save filtered data
% parfor i=1:S1(4)
% Dnewtest(:,:,:,i)=reshape(Rnewtest(i,:),[S1(1),S1(2),S1(3)]);
% end
% out4.vol = Dnewtest;
% err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCAFtest.nii'));




