clear all

% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);

load(strcat(file(1:size(file,2)-16),'TS.mat'))

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

imagesc(Amod*10)

TSav=zeros(max(M),size(TS,2));   

for i=1:max(M)
    Dm=D(:,:,:,M==i);
    out.vol = Dm;
    err = MRIwrite(out,strcat(file(1:size(file,2)-6),'diffmodule',int2str(i),'.nii'));
    DlM(:,:,:,i)=sum(Dm,4);
    TSm=TS(M==i,:);
    TSav(i,:)=mean(TS(M==i,:));
    save(strcat(file(1:size(file,2)-10),'TSmodule',num2str(i),'.mat'),'TSm');
end

save(strcat(file(1:size(file,2)-6),'TSav.mat'),'TSav');

DlM(DlM>1)=1;   
Dcolor=zeros(S1(1),S1(2),S1(3),3);
C(:,1)=[1,0,0];
C(:,2)=[0,1,0];
C(:,3)=[0,0,1];
C(:,4)=[1,1,0];
C(:,5)=[0,1,1];
C(:,6)=[1,0,1];
 
for i=1:max(M)
    Dcolor(:,:,:,1)=Dcolor(:,:,:,1)+DlM(:,:,:,i)*C(1,(i-ceil((i-6)/6)*6));
    Dcolor(:,:,:,2)=Dcolor(:,:,:,2)+DlM(:,:,:,i)*C(2,(i-ceil((i-6)/6)*6));
    Dcolor(:,:,:,3)=Dcolor(:,:,:,3)+DlM(:,:,:,i)*C(3,(i-ceil((i-6)/6)*6));
end

Dcolor(:,:,:,1)=Dcolor(:,:,:,1)/(max(max(max(Dcolor(:,:,:,1)))));
Dcolor(:,:,:,2)=Dcolor(:,:,:,2)/(max(max(max(Dcolor(:,:,:,2)))));
Dcolor(:,:,:,3)=Dcolor(:,:,:,3)/(max(max(max(Dcolor(:,:,:,3)))));

D1m=Montage3(Dcolor(:,:,:,1));
D2m=Montage3(Dcolor(:,:,:,2));
D3m=Montage3(Dcolor(:,:,:,3));
D2m=D2m/max(max(D2m));
D1m=D1m/max(max(D1m));
D3m=D3m/max(max(D3m));



Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); 
Pc=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'coloredSumMapdiff.png'));
imwrite(Dm/Pc,fullFileName);