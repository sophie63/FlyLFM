
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


% Demean
S2=size(R);
parfor i=1:S2(2)
    Rm(i)=mean(R(:,i));
    R(:,i)=R(:,i)-Rm(i);
end

P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
R1=R.*(10000/P);

R=R1;

% % %simple std norm (choose between this and Smith lines below)
%   stddevs=max(std(R),0.01);  
%   R=R./repmat(stddevs,size(R,1),1);  % var-norm

%%Variance normalisation ala melodic
% [uu,ss,vv]=nets_svds(R,30); % initial SVD to the top 30 components (arbitrary number fixed in MELODIC)
% %vv(abs(vv)<std(vv(:)))=0;
% stddevs=max(std(R-uu*ss*vv'),1);  % subtract main parts of top components from data to get normalisation
% R=R./repmat(stddevs,size(R,1),1);  % var-norm

%SVD
Npc=60;
[u,s,v]=nets_svds(R,Npc);


%% ICA
[icasig, A, W] = fastica (v','approach','symm','epsilon', 0.001);

% form back ica maps
GM=v*A;
GMv=GM./repmat(std(GM),size(GM,1),1);
%Correct sign so that mean of the positive side is larger than the mean of
%the negative side after zscoring

GMzp=GMv-2;
GMzp(GMzp<0)=0;
GMzpm=mean(GMzp);
GMzn=GMv+2;
GMzn(GMzn>0)=0;
GMznm=mean(GMzn);
GMs=GM.*repmat(sign(GMzpm+GMznm),size(GM,1),1); 

TS=u*s*A;
TSs=TS.*repmat(sign(GMzpm+GMznm),size(TS,1),1);


% Reform volumes
for i=1:Npc
Dica(:,:,:,i)=reshape(GMs(:,i),[S1(1),S1(2),S1(3)]);
end

% Save ICA maps and time series
out.vol = Dica;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'IC.nii'));

save(strcat(file(1:size(file,2)-4),'TS'),'TSs')

pause

% Now open in Sort_ICs_from_matlab_QUICK_fortemp.ipynb to find activity
% components

% %% Reform from non activity components
% Rfromica=u(:,I)*s(I,I)*A(I,I)*icasig(I,:);
% Rfromicab=Rfromica.*repmat(stddevs,size(R,1),1);
% clear Rfromica
% R
% Rfromica2=Rfromicab.*P/10000;
% clear Rfromicab
% parfor i=1:S2(2)
%     Rb(:,i)=Rfromica2(:,i)+Rm(i);
% end
% clear Rfromica2
 
%% Reform removing activity components
Rfromica=u*s*A(:,I)*icasig(I,:);
 parfor i=1:S1(4)
 NewDfromicaa(:,:,:,i)=reshape(Rfromica(i,:),[S1(1),S1(2),S1(3)]);
 end
out4.vol = NewDfromicaa;
err = MRIwrite(out4,'/media/Buffer3/862shortonlyactivity.nii');

Rc=R-Rfromica;
 parfor i=1:S2(2)
    Rb(:,i)=Rc(:,i)+Rm(i).*(10000/P);
 end
 clear Rfromica
 
%% Save template(R(:,i)-Rfromica(:,i)).*repmat(stddevs,size(R,1),1)
 parfor i=1:S1(4)
 NewDfromica(:,:,:,i)=reshape(Rb(i,:),[S1(1),S1(2),S1(3)]);
 end
 
 clear Rb
 
 out5.vol = NewDfromica;
 err = MRIwrite(out5,'/media/Buffer3/862shorttemplate.nii');
