clear all

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);

load(strcat(file(1:size(file,2)-6),'TS.mat'))

TS=TSo';
Acov=TS*TS';
 Ainv=pinv(Acov);
DTS=diff(TS');
 DTS2=[DTS' DTS(end,:)'];
 Apart=DTS2*TS';
 Afinal=Apart*Ainv;
 Adiff=DTS2*DTS2';

W=Adiff;
n  = size(W,1);             % number of nodes
M  = 1:n;                   % initial community affiliations
Q0 = -1; Q1 = 0;            % initialize modularity values
while Q1-Q0>1e-5;           % while modularity increases
Q0 = Q1;                % perform community detection
%[M, Q1] = community_louvain(W, [], M,'negative_sym');
[M, Q1] = community_louvain(W, [], M,'negative_asym');
end
imagesc(Adiff)
[Or Amod]=reorder_mod(W,M);

imagesc(Amod)

for i=1:max(M)
    Dm=D(:,:,:,M==i);
    size(Dm)
    if size(size(Dm),2)==4
        ColoredSumaryMapfunc(strcat(file(1:size(file,2)-6),'module',i,'.nii'),Dm);
end
end
    

%ModuleIndices=M(