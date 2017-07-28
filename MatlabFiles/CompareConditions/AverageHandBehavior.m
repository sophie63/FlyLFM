%choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

%% Average data during walk and during rest
Dwalk=zeros(S(1),S(2),S(3),1);
Drest=zeros(S(1),S(2),S(3),1);
Drgroom=zeros(S(1),S(2),S(3),1);
Dfgroom=zeros(S(1),S(2),S(3),1);
Dgroom=zeros(S(1),S(2),S(3),1);

for i =1:S(4)
    Dwalk=Dwalk+Data(:,:,:,i)*WalkHandk(i);
    Drest=Drest+Data(:,:,:,i)*RestHandk(i);
    %Dfgroom=Dfgroom+Data(:,:,:,i)*FrontGroomk(i);
    %Drgroom=Drgroom+Data(:,:,:,i)*RearGroomk(i);
    Dgroom=Dgroom+Data(:,:,:,i)*Groomk(i);
end

Dwalk=Dwalk/sum(WalkHandk);
Drest=Drest/sum(RestHandk);
%Dfgroom=Dfgroom/sum(FrontGroomk);
%Drgroom=Drgroom/sum(RearGroomk);
Dgroom=Dgroom/sum(Groomk);

DrestN=Drest-Dwalk;
DrestN(DrestN<0)=0;
DwalkN=Dwalk-Drest;
DwalkN(DwalkN<0)=0;

out.vol=cat(4,DrestN,DwalkN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestWalk.nii'));
 
% DrestN=Drest-Dfgroom;
% DrestN(DrestN<0)=0;
% DfgroomN=Dfgroom-Drest;
% DfgroomN(DfgroomN<0)=0;
% 
% out.vol=cat(4,DrestN,DfgroomN);
% err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestGroomFront.nii'));
% 
% DrestN=Drest-Drgroom;
% DrestN(DrestN<0)=0;
% DrgroomN=Drgroom-Drest;
% DrgroomN(DrgroomN<0)=0;
% 
% out.vol=cat(4,DrestN,DrgroomN);
% err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestGroomBack.nii'));

DrestN=Drest-Dwalk;
DrestN(DrestN<0)=0;
DwalkN=Dwalk-Drest;
DwalkN(DwalkN<0)=0;

out.vol=cat(4,DrestN,DwalkN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestWalk.nii'));

DrestN=Drest-Dgroom;
DrestN(DrestN<0)=0;
DgroomN=Dgroom-Drest;
DgroomN(DgroomN<0)=0;

out.vol=cat(4,DrestN,DgroomN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestGroom.nii'));

DgroomN=Dgroom-Dwalk;
DgroomN(DgroomN<0)=0;
DwalkN=Dwalk-Dgroom;
DwalkN(DwalkN<0)=0;

out.vol=cat(4,DgroomN,DwalkN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'GroomWalk.nii'));

