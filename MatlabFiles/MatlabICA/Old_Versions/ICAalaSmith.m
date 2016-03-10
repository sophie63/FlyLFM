clear all

%open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D1=B.vol;

% %grand mean intensity normalisation
% Dm=squeeze(min(min(min(min(D1)))));
% D2=D1-Dm;
% DM=squeeze(max(max(max(max(D2)))));
% D=D2./DM;

D=D1;

%Masking of non brain data
%get the data before unbleaching

%reshape
S1=size(D);

% %boxav
% for i=1:S1(3)
%     for j=1:S1(4)
%         X=squeeze(D(1:2*(floor(S1(1)/2)),1:2*(floor(S1(2)/2)),i,j));
%         Xs=reshape(sum(sum(reshape(X, 2, floor(S1(1)/2), 2, floor(S1(2)/2)), 1), 3), floor(S1(1)/2), floor(S1(2)/2));
%         Ds(:,:,i,j)=Xs;
%     end
%     i
% end
% 
% D=Ds;

S1=size(D);
% 
% for i=1:S1(3)
% C=squeeze(D(:,:,i,:));
% k75=Kalman_Stack_Filter(C,0.75);
% Dkf(:,:,i,:)=k75;
% end
% 
% D=Dkf;

parfor i=1:S1(4)
R(:,i)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end

%add floats or double here?
P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=double(R).*(10000/P);

% demean
S2=size(R1);
parfor i=1:S2(1)
    R2(i,:)=R1(i,:)-mean(R1(i,:));
    %R2(i,:)=medfilt3(R2(i,:),5);
end

%need timeXspace
R=R2';

% % %simple std norm (choose between this and Smith lines below)
  stddevs=max(std(R),0.001);  
  R=R./repmat(stddevs,size(R,1),1);  % var-norm

%%Smith lines
%[uu,ss,vv]=nets_svds(R,30); % initial SVD to the top 30 components (arbitrary number fixed in MELODIC)
%vv(abs(vv)<2.3*std(vv(:)))=0;  % zero small values in PCA maps  %what is that doing here?
%stddevs=max(std(R-uu*ss*vv'),0.1);  % subtract main parts of top components from data to get normalisation
%R=R./repmat(stddevs,size(R,1),1);  % var-norm

%SVD
Npc=30;
[u,s,v]=nets_svds(R,Npc);

% can be replaced by matlab svd : 
%[u,s,v]=svd(R);
%u=u(:,1:Npc);
%v=v(:,1:Npc);
%s=s(1:Npc,1:Npc);

%5/9/2015: add z scoring before ICA
%v(abs(v)<2.3*std(v(:)))=0; 
mixedsig = s * v';


%pseudo time series for fsl
DE=s*v';
for i=1:Npc
Dpca(:,:,:,i)=reshape(DE(i,:),[S1(1),S1(2),S1(3)]);
end

% save PCA maps
out4.vol = Dpca;
err = MRIwrite(out4,'/home/sophie/Desktop/Av863_sv40.nii');


%ICA - using v instead of mixedsign == whitening
[icasig, A, W] = fastica (v','approach','symm','epsilon', 0.000001);
%check again difference whithout symm

% in fsl icaobj.perf_ica(melodat.get_white()*melodat.get_Data()); v*s*v'?

% form back ica maps
GM=v*A;

% note that *here* GM is spaceXcomponents
GM=GM.*repmat(sign(max(GM)+min(GM)),size(GM,1),1)./repmat(std(GM),size(GM,1),1); 
% make all individual group maps have a positive peak, and of peak height=1

%TS=u(melodicmix;
TS=TS.*repmat(sign(max(GM)+min(GM)),size(TS,1),1);

for i=1:Npc
Dica(:,:,:,i)=reshape(GM(:,i),[S1(1),S1(2),S1(3)]);
end

figure
show_IC(Dica)

save ICA maps
out.vol = Dica;
err = MRIwrite(out,'C:\Users\admin\Desktop\862ss2c500regcUKF_ICs30.nii');

% for i=1:Npc
% svdresh(:,:,:,i)=reshape(v(:,i),[S1(1),S1(2),S1(3)]);
% end
% out.vol = svdresh;
% err = MRIwrite(out,'C:\Users\admin\Desktop\vresh.nii');

