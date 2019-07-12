clear all

% Get regressors
[FileName,PathName] = uigetfile('*.mat','Select the Xk file','/media/sophie2/');
file=strcat(PathName,FileName)
load(file);
X=Xk;
Sx=size(X);

% Choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);


%% Average data during walk, groom, rest and turn
Dpre=zeros(S(1),S(2),S(3),1);
Dpanic=zeros(S(1),S(2),S(3),1);
Drelief=zeros(S(1),S(2),S(3),1);

Dwalk=zeros(S(1),S(2),S(3),1);
for i =1:S(4)
    Dpre(:,:,:,1)=Dpre(:,:,:,1)+Data(:,:,:,i)*X(1,i);
    Dpanic(:,:,:,1)=Dpanic(:,:,:,1)+Data(:,:,:,i)*X(2,i);    
    Drelief(:,:,:,1)=Drelief(:,:,:,1)+Data(:,:,:,i)*X(3,i);
    Dwalk(:,:,:,1)=Dwalk(:,:,:,1)+Data(:,:,:,i)*X(4,i);
end

Dpre=Dpre/sum(X(1,:));
Dpanic=Dpanic/sum(X(2,:));
Drelief=Drelief/sum(X(3,:));
Dwalk=Dwalk/sum(X(4,:));

DpanicN=Dpanic-Dpre;
D2=DpanicN;
D2(DpanicN<0)=0;
D1=-DpanicN;
D1(DpanicN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'PreFlail.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'PreFlailM.PNG'));
imwrite(Dm/M, fullFileName);



DreliefN=Drelief-Dpre;
D2=DreliefN;
D2(DreliefN<0)=0;
D1=-DreliefN;
D1(DreliefN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'PreRelief.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'PreReliefM.PNG'));
imwrite(Dm/M, fullFileName);

DpanicN=Dpanic-Drelief;
D2=DpanicN;
D2(DpanicN<0)=0;
D1=-DpanicN;
D1(DpanicN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'ReliefFlail.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'ReliefFlailM.PNG'));
imwrite(Dm/M, fullFileName);


DpanicN=Dpanic-Dwalk;
D2=DpanicN;
D2(DpanicN<0)=0;
D1=-DpanicN;
D1(DpanicN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'WalkFlail.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'WalkFlailM.PNG'));
imwrite(Dm/M, fullFileName);




