clear all

% Get regressors
[FileName,PathName] = uigetfile('*.mat','Select the Xk file','/media/sophie2/');
file=strcat(PathName,FileName)
load(file);
Sx=size(Xk);

% Choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

% Plot and choose a threshold
plot(Xk(:,4))

Thresh=200;

% Rest: if not walk, not groom (and not weird if that would be there)
Walkt=zeros(Sx(1),1);
Walkt(Xk(:,3)>Thresh)=1;

Groomt=zeros(Sx(1),1);
Groomt((Xk(:,5)>50)&(Xk(:,3)<Thresh))=1;

Rest=zeros(Sx(1),1);
Rest(~(logical(Groomt)|logical(Walkt)))=1;

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
    Dleft(:,:,:,1)=Dleft(:,:,:,1)+Data(:,:,:,i)*Xk(i,1);
    Dright(:,:,:,1)=Dright(:,:,:,1)+Data(:,:,:,i)*Xk(i,2);
end

Dwalk=Dwalk/sum(Walkt);
Dgroom=Dgroom/sum(Groomt);
Drest=Drest/sum(Rest);
Dleft=Dleft/sum(Xk(:,1));
Dright=Dright/sum(Xk(:,2));

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
M=max(max(max(max(Dm4norm))));
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
M=max(max(max(max(Dm4norm))));
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
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'RestGroomM.PNG'));
imwrite(Dm/M, fullFileName);

Dmgroom=Dm;

DmGW=cat(2,DmWalk,Dmgroom);
Dm4norm=DmGW;
Dm4norm(DmGW==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'WalkGroomM.PNG'));
imwrite(DmGW/M, fullFileName);


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
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'GroomWalkM.PNG'));
imwrite(Dm/M, fullFileName);


%% Now measure average and variance

Dwalk4v=Data(:,:,:,Walkt==1);
Dgroom4v=Data(:,:,:,Groomt==1);
Drest4v=Data(:,:,:,Rest==1);

AvWalk=mean(reshape(Dwalk4v,S(1)*S(2)*S(3),size(Dwalk4v,4)),1);
AvAvWalk=mean(AvWalk);
VarAvWalk=var(AvWalk);

AvGroom=mean(reshape(Dgroom4v,S(1)*S(2)*S(3),size(Dgroom4v,4)),1);
AvAvGroom=mean(AvGroom);
VarAvGroom=var(AvGroom);

AvRest=mean(reshape(Drest4v,S(1)*S(2)*S(3),size(Drest4v,4)),1);
AvAvRest=mean(AvRest);
VarAvRest=var(AvRest);

fullFileName = fullfile(strcat(file(1:size(file,2)-4),'AvVar.mat'));
save(fullFileName,'AvAvWalk','VarAvWalk','AvAvGroom','VarAvGroom','AvAvRest','VarAvRest')

fullFileName = fullfile(strcat(file(1:size(file,2)-4),'AvVar.txt'));
fileID = fopen(fullFileName,'w');
fprintf(fileID,fullFileName,'\r\n');
fprintf(fileID,'AvAvWalk,VarAvWalk,AvAvGroom,VarAvGroom,AvAvRest,VarAvRest\r\n');
fprintf(fileID,'%12.8f\r\n',AvAvWalk);
fprintf(fileID,'%12.8f\n',VarAvWalk);
fprintf(fileID,'%12.8f\n',AvAvGroom);
fprintf(fileID,'%12.8f\n',VarAvGroom);
fprintf(fileID,'%12.8f\n',AvAvRest);
fprintf(fileID,'%12.8f\n',VarAvRest);
fclose(fileID);

