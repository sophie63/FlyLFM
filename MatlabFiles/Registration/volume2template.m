clear all
tic

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=B.vol;

% Reshape
S1=size(D);
parfor i=1:S1(4)
R(:,i)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end

% demean
S2=size(R);
parfor i=1:S2(1)
    R2(i,:)=R(i,:)-mean(R(i,:));
end

% Need timeXspace
R=R2';

% SVD
Npc=30;
[u,s,v]=nets_svds(R,Npc);

% Save PCA components to look at in FIJI
DE=s*v';
for i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end

% save PCA maps
out4.vol = Dpca;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCA.nii'));
toc

pause

% Remove activity component
NewR=u*s*v;

% Remean
parfor i=1:S2(1)
    NewR2(i,:)=NewR(i,:)+mean(R1(i,:));
end

% Save volume template
parfor i=1:S1(4)
    NewD2(:,:,:,i)=reshape(NewR2(i,:),[S1(1),S1(2),S1(3)]);
end
out.vol = NewD2;
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'temp.nii'));

