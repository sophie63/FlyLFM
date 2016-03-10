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

D=double(D1);

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
R1=R.*(10000/P);

% demean
S2=size(R1);
parfor i=1:S2(1)
    R2(i,:)=R1(i,:)-mean(R1(i,:));
    %R2(i,:)=medfilt3(R2(i,:),5);
end

%need timeXspace
R=R2';

% % %simple std norm (choose between this and Smith lines below)
%  stddevs=max(std(R),0.01);  
%  R=R./repmat(stddevs,size(R,1),1);  % var-norm

%%Smith lines
[uu,ss,vv]=nets_svds(R,30); % initial SVD to the top 30 components (arbitrary number fixed in MELODIC)
%vv(abs(vv)<2.3*std(vv(:)))=0;  % zero small values in PCA maps  
vv(abs(vv)<0.23*std(vv(:)))=0;
stddevs=max(std(R-uu*ss*vv'),1);  % subtract main parts of top components from data to get normalisation
R=R./repmat(stddevs,size(R,1),1);  % var-norm

%SVD
Npc=200;
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
err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCA.nii'));

%ICA - using v instead of mixedsign == whitening
%[icasig, A, W] = fastica (v','approach','symm','epsilon', 0.000000001);
[icasig, A, W] = fastica (v','epsilon', 0.01);

% in fsl icaobj.perf_ica(melodat.get_white()*melodat.get_Data()); v*s*v'?

% form back ica maps
GM=v*A;

%correct sign
vartry=sqrt(var(GM));
for i=1:S2(1)
GMv(i,:)=GM(i,:)./vartry;
end
GMzp=GMv;
GMzp(GMzp<3)=0;
GMzpm=mean(GMzp);
GMzn=GMv;
GMzn(GMzn>-3)=0;
GMznm=mean(GMzn);
% note that *here* GM is spGMznmedaceXcomponents
GMs=GM.*repmat(sign(GMzpm+GMznm),size(GM,1),1)./repmat(std(GM),size(GM,1),1); 
% make all individual group maps have a positive peak, and of peak height=1
%could use sign(max(median(zscoredpos),(zscoredneg)) instead

TS=u*s*A;
TSs=TS.*repmat(sign(GMzpm+GMznm),size(TS,1),1);


%reorder maps by variance
[varo,Order]=sort(var(TS));
TSo=TSs(:,Order);
GMo=GMs(:,Order);

for i=1:Npc
Dica(:,:,:,i)=reshape(GMs(:,i),[S1(1),S1(2),S1(3)]);
end
% in fsl icaobj.perf_ica(melodat.get_white()*melodat.get_Data()); v*s*v'?


%err = MRIwrite(out4,strcat(file(1:size(file,2)-4),'PCA.nii')save ICA maps
out.vol = Dica;
%out.vol=D*100000/P;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'IC.nii'))

save(strcat(file(1:size(file,2)-4),'TS'),'TSs')

% %save zscored ica maps
% Dicaz=Dica;
% Dicaz(Dicaz<2)=2;
% Dicaz=Dicaz-2;
% Dicabin=Dicaz;
% Dicabin(Dicabin>0)=1;
% 
% 
% %Calculate antimap to correct for everything below z scored threshold
% Dicaantiz=GM;
% Dicaantiz(Dicaantiz>2)=0;
% 
% Correction=zeros(S2(2),Npc);
% 
% 
% % Reconstruct the signal
% NewR=TS*GM';
% parfor i=1:S1(4)
% NewD(:,:,:,i)=reshape(NewR(i,:),[S1(1),S1(2),S1(3)]);
% end
% 
% out3.vol = NewD;
% %out.vol=D*100000/P;
% err = MRIwrite(out3,'/home/sophie/Desktop/863rec.nii');
% 
% %Only component 2
% NewR2=TS(:,2)*GM(:,2)';
% parfor i=1:S1(4)
% NewD2(:,:,:,i)=reshape(NewR2(i,:),[S1(1),S1(2),S1(3)]);
% end
% out4.vol = NewD2;
% %out.vol=D*100000/P;
% err = MRIwrite(out4,'/home/sophie/Desktop/863reccomp2.nii');
% 
% 
% Rfromica=u*s*A*icasig;
% parfor i=1:S1(4)
% NewDfromica(:,:,:,i)=reshape(Rfromica(i,:),[S1(1),S1(2),S1(3)]);
% end
% out5.vol = NewDfromica;
% %out.vol=D*100000/P;
% err = MRIwrite(out5,'/home/sophie/Desktop/863fromica.nii');