<<<<<<< HEAD
% This script performs PCA and spatial ICA in a manner similar to melodic from FSL 

clear all

% Open Data
InputDir= uigetdir;
files=dir(InputDir)

for k=3:length(files)
B=MRIread(strcat(InputDir,'/',files(k).name));
D=double(B.vol);
clear B

% Reshape as a space*time 2D matrix 
S1=size(D);
parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end
clear D


% Demean
S2=size(R);
parfor i=1:S2(2)
    Rm(i)=mean(R(:,i));
    R(:,i)=R(:,i)-Rm(i);
end

P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=R.*(10000/P);
clear R
R=R1;

% % %simple std norm (choose between this and Smith lines below)
   stddevs=max(std(R),0.01);  
   R=R./repmat(stddevs,size(R,1),1);  % var-norm



%% Variance normalisation ala melodic
%[uu,ss,vv]=nets_svds(R,10); 
% initial SVD to the top components
%vv(abs(vv)<0.4*std(vv(:)))=0;
%vv(abs(vv)<2.3*std(vv(:)))=0;
%vv(abs(vv)<0.023*std(vv(:)))=0;
%stddevs=max(std(R-uu*ss*vv'),1);  
% Subtract main parts of top components from data to get normalisation
% Which leaves most main components not normalized
% Note also that small values are not increased (normalized by one)
%R=R./repmat(stddevs,size(R,1),1);  % Var-norm

% Check dimentionality
[uu,ss,vv]=nets_svds(R,3000); 
plot(log(diag(ss)))

file=strcat(InputDir,'/',files(k).name)

file2=strcat(file(1:size(file,2)-4),'Dim.mat')
var=diag(ss);
save(file2,'var')

end
||||||| merged common ancestors
=======
% This script performs PCA and spatial ICA in a manner similar to melodic from FSL 

clear all

% Open Data
InputDir= uigetdir;
files=dir(InputDir)

for k=3:length(files)
B=MRIread(strcat(InputDir,'/',files(k).name));
D=double(B.vol);
clear B

% Reshape as a space*time 2D matrix 
S1=size(D);
parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end
clear D


% Demean
S2=size(R);
parfor i=1:S2(2)
    Rm(i)=mean(R(:,i));
    R(:,i)=R(:,i)-Rm(i);
end

P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=R.*(10000/P);
clear R
R=R1;

% % %simple std norm (choose between this and Smith lines below)
   stddevs=max(std(R),0.01);  
   R=R./repmat(stddevs,size(R,1),1);  % var-norm



%% Variance normalisation ala melodic
%[uu,ss,vv]=nets_svds(R,10); 
% initial SVD to the top components
%vv(abs(vv)<0.4*std(vv(:)))=0;
%vv(abs(vv)<2.3*std(vv(:)))=0;
%vv(abs(vv)<0.023*std(vv(:)))=0;
%stddevs=max(std(R-uu*ss*vv'),1);  
% Subtract main parts of top components from data to get normalisation
% Which leaves most main components not normalized
% Note also that small values are not increased (normalized by one)
%R=R./repmat(stddevs,size(R,1),1);  % Var-norm

% Check dimentionality
[uu,ss,vv]=nets_svds(R,3000); 
plot(log(diag(ss)))

file=strcat(InputDir,'/',files(k).name)

file2=strcat(file(1:size(file,2)-4),'Dim.mat')
var=diag(ss);
save(file2,'var')
clear D
clear S1
clear R
clear R1
clear uu
clear ss
clear vv
end
>>>>>>> 57759d29a73933bc92c6aaedf808d46b8fabe631
