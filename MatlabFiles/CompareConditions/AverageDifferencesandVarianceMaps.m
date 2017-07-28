clear all

% Get regressors
[FileName,PathName] = uigetfile('*.mat','Select the Xk file','/media/sophie2/');
file=strcat(PathName,FileName)
load(file);
%Xk=Xk';
Sx=size(Xk);


% Choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

% Plot and choose a threshold
plot(Xk(:,3))
prompt = 'What is the threshold for walking?';
Thresh = input(prompt)
%Thresh=0.2;

% Rest: if not walk, not groom (and not weird if that would be there)
Walkt=zeros(Sx(1),1);
Walkt(Xk(:,3)>Thresh)=1;

Groomt=zeros(Sx(1),1);
%Groomt((Xk(:,5)>200)&(Xk(:,4)<Thresh))=1;

Rest=zeros(Sx(1),1);
Rest(~(logical(Groomt)|logical(Walkt)))=1;

LeftminRight=Xk(:,1)-Xk(:,2);
Left2=LeftminRight;
Left2(LeftminRight<0)=0;
Right2=-LeftminRight;
Right2(LeftminRight>0)=0;

%% Average data during walk, groom, rest and turn
Dwalk=zeros(S(1),S(2),S(3),1);
Dgroom=zeros(S(1),S(2),S(3),1);
Drest=zeros(S(1),S(2),S(3),1);
Dleft=zeros(S(1),S(2),S(3),1);
Dright=zeros(S(1),S(2),S(3),1);

for i =1:S(4)
    Dwalk(:,:,:,1)=Dwalk(:,:,:,1)+Data(:,:,:,i)*Walkt(i);
    Dgroom(:,:,:,1)=Dgroom(:,:,:,1)+Data(:,:,:,i)*Groomt(i);    
    Drest(:,:,:,1)=Drest(:,:,:,1)+Data(:,:,:,i)*Rest(i);
    Dleft(:,:,:,1)=Dleft(:,:,:,1)+Data(:,:,:,i)*Left2(i);
    Dright(:,:,:,1)=Dright(:,:,:,1)+Data(:,:,:,i)*Right2(i);
end

Dwalk=Dwalk/sum(Walkt);
Dgroom=Dgroom/sum(Groomt);
Drest=Drest/sum(Rest);
Dleft=Dleft/sum(Left2);
Dright=Dright/sum(Right2);

DleftrightN=Dleft-Dright;
D1=DleftrightN;
D1(DleftrightN<0)=0;
D2=-DleftrightN;
D2(DleftrightN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'LeftMinRight.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'LeftMinRightM.PNG'));
imwrite(Dm/M, fullFileName);

DWalkN=Dwalk-Drest;
D2=DWalkN;
D2(DWalkN<0)=0;
D1=-DWalkN;
D1(DWalkN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestWalk.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'RestWalkM.PNG'));
imwrite(Dm/M, fullFileName);

DmWalk=Dm;

DGroomN=Dgroom-Drest;
D2=DGroomN;
D2(DGroomN<0)=0;
D1=-DGroomN;
D1(DGroomN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestGroom.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'RestGroomM.PNG'));
imwrite(Dm/M, fullFileName);

Dmgroom=Dm;

DmGW=cat(2,DmWalk,Dmgroom);
Dm4norm=DmGW;
Dm4norm(DmGW==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'WalkGroomM.PNG'));
imwrite(DmGW*2/M, fullFileName);


DWalkN=Dwalk-Dgroom;
D2=DWalkN;
D2(DWalkN<0)=0;
D1=-DWalkN;
D1(DWalkN>0)=0;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'GroomWalk.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'GroomWalkM.PNG'));
imwrite(Dm/M, fullFileName);


%% Now make maps of variance

Dwalk4v=Data(:,:,:,Walkt==1);
Dgroom4v=Data(:,:,:,Groomt==1);
Drest4v=Data(:,:,:,Rest==1);
Dleftfull=Data(:,:,:,Left2>6*10^4);
Drightfull=Data(:,:,:,Right2>6*10^4);

out.vol=Dleftfull;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'OnlyLeft.nii'));

out.vol=Drightfull;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'OnlyRight.nii'));

out.vol=Dwalk4v;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'OnlyWalk.nii'));

out.vol=Dgroom4v;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'OnlyGroom.nii'));

out.vol=Drest4v;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'OnlyRest.nii'));

Dvarwalk=zeros(S(1),S(2),S(3),1);
Dvargroom=zeros(S(1),S(2),S(3),1);
Dvarrest=zeros(S(1),S(2),S(3),1);

for i=1:S(1)
    for k=1:S(3)
       parfor j=1:S(2)        
Dvarwalk(i,j,k)=var(Dwalk4v(i,j,k,:));
Dvargroom(i,j,k)=var(Dgroom4v(i,j,k,:));
Dvarrest(i,j,k)=var(Drest4v(i,j,k,:));
        end
    end
    i
end



out.vol=cat(4,Dvarrest,Dvarwalk);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestWalkVar.nii'));

D1=Dvarrest;
D2=Dvarwalk;
D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'RestWalkVarM.PNG'));
imwrite(Dm/M, fullFileName);


D1=Dvarrest;
D2=Dvargroom;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestGroomVar.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'RestGroomVarM.PNG'));
imwrite(Dm/M, fullFileName);
