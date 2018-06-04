%choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

%Average walk and not walk
%Walk=Leftk+Rightk+Straightk;
Walk=Xk(4:7249,3);
Rest=Xk(4:7249,6);
Groom=Xk(4:7249,5);
%plot and choose a threshold
plot(Walk)

%% Average data during walk and during rest
Dwalk=zeros(S(1),S(2),S(3),1);
Drest=zeros(S(1),S(2),S(3),1);
Dstraight=zeros(S(1),S(2),S(3),1);
Dturn=zeros(S(1),S(2),S(3),1);
Dleft=zeros(S(1),S(2),S(3),1);
Dright=zeros(S(1),S(2),S(3),1);
Dgroom=zeros(S(1),S(2),S(3),1);
L=0;
R=0;

for i =1:S(4)
    Dwalk=Dwalk+Data(:,:,:,i)*Walk(i);
    Drest=Drest+Data(:,:,:,i)*Rest(i);
    Dgroom=Dgroom+Data(:,:,:,i)*Groom(i);
%    Dstraight=Dstraight+Data(:,:,:,i)*Straightk(i);
%    Dturn=Dturn+Data(:,:,:,i)*(Leftk(i)+Rightk(i));
%    if (Leftk(i)>Rightk(i))
%    Dleft=Dleft+Data(:,:,:,i)*(Leftk(i)-Rightk(i));
%    L=L+Leftk(i)-Rightk(i);
%    elseif (Rightk(i)>Leftk(i))
%    Dright=Dright+Data(:,:,:,i)*(Rightk(i)-Leftk(i));
 %   R=R+(Rightk(i)-Leftk(i));
end

Dwalk=Dwalk/sum(Walk);
Drest=Drest/sum(Rest);

%Dleft=Dleft/L;
%Dright=Dright/R;
%Dstraight=Dstraight/sum(Straightk);
%Dturn=Dturn/sum(Leftk+Rightk);

DrestN=Drest-Dwalk;
DrestN(DrestN<0)=0;
DwalkN=Dwalk-Drest;
DwalkN(DwalkN<0)=0;

out.vol=cat(4,DrestN,DwalkN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestWalk.nii'));

Dgroom=Dgroom/sum(Groom);
Drest=Drest/sum(Rest);
DrestN=Drest-Dgroom;
DrestN(DrestN<0)=0;
DwalkN=Dwalk-Drest;
DwalkN(DwalkN<0)=0;

out.vol=cat(4,DrestN,DwalkN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestGroom.nii'));


% DstraightN=Dstraight-Dturn;
% DstraightN(DstraightN<0)=0;
% DturnN=Dturn-Dstraight;
% DturnN(DturnN<0)=0;
% out.vol=cat(4,DstraightN,DturnN);
% err = MRIwrite(out,strcat(file(1:size(file,2)-4),'StraightTurn.nii'));
% 
% DleftN=Dleft-Dright;
% DleftN(DleftN<0)=0;
% DrightN=Dright-Dleft;
% DrightN(DrightN<0)=0;
% out.vol=cat(4,DleftN,DrightN);
% err = MRIwrite(out,strcat(file(1:size(file,2)-4),'LeftRight.nii'));
