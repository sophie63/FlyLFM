%open oxygen 100081 aligned

%choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

DataDown=zeros(S(1),S(2),S(3),2800);
DataUp=zeros(S(1),S(2),S(3),2800);
DataTrans=zeros(S(1),S(2),S(3),240);
DataDownStrong=zeros(S(1),S(2),S(3),1400);
DataUpStrong=zeros(S(1),S(2),S(3),1400);
DataStrongHyp=zeros(S(1),S(2),S(3),1400);

for i=1:4
    DataDown(:,:,:,((i-1)*700+1:i*700))=Data(:,:,:,(i*1800+101):(i*1800+900)-100);
    DataUp(:,:,:,((i-1)*700+1:i*700))=Data(:,:,:,(i*1800+101+900):(i*1800+1800)-100);
    DataTrans(:,:,:,((i-1)*60+1:i*60))=Data(:,:,:,(i*1800-29):(i*1800+30)); 
    DataDownStrong(:,:,:,((i-1)*350+1:i*350))=Data(:,:,:,(i*1800+100):(i*1800+449));
    DataUpStrong(:,:,:,((i-1)*350+1:i*350))=Data(:,:,:,(i*1800+100+900):(i*1800+900+449));
    DataStrongHyp(:,:,:,((i-1)*350+1:i*350))=Data(:,:,:,(i*1800+450):(i*1800+799));    
end

DDown=sum(DataDown,4)/size(DataDown,4);
DUp=sum(DataUp,4)/size(DataUp,4);
DDownStrong=sum(DataDownStrong,4)/size(DataDownStrong,4);
DUpStrong=sum(DataUpStrong,4)/size(DataUpStrong,4);
DTrans=sum(DataTrans,4)/size(DataTrans,4);
DStrongHyp=sum(DataStrongHyp,4)/size(DataStrongHyp,4);

DataControl=Data(:,:,:,300:1600);
DataHyp=Data(:,:,:,1900:end);

DControl=sum(DataControl,4)/size(DataControl,4);
DHyp=sum(DataHyp,4)/size(DataHyp,4);

DControlN=DControl-DHyp;
DControlN(DControlN<0)=0;
DHypN=DHyp-DControl;
DHypN(DHypN<0)=0;

out.vol=cat(4,DControlN,DHypN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'ControlHypN.nii'));

DDownN=DDown-DUp;
DDownN(DDownN<0)=0;
DUpN=DUp-DDown;
DUpN(DUpN<0)=0;

out.vol=cat(4,DDownN,DUpN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'DownUpN.nii'));

DDownStrongN=DDownStrong-DUpStrong;
DDownStrongN(DDownStrongN<0)=0;
DUpStrongN=DUpStrong-DDownStrong;
DUpStrongN(DUpStrongN<0)=0;

out.vol=cat(4,DDownStrongN,DUpStrongN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'DownUpStrongN.nii'));

DControlN=DControl-DTrans;
DControlN(DControlN<0)=0;
DTransN=DTrans-DControl;
DTransN(DTransN<0)=0;

out.vol=cat(4,DControlN,DTransN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'ControlTransitionN.nii'));

DControlN=DControl-DStrongHyp;
DControlN(DControlN<0)=0;
DStrongHypN=DStrongHyp-DControl;
DStrongHypN(DStrongHypN<0)=0;


out.vol=cat(4,DControlN,DStrongHypN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'ControlStrongHypoxiaN.nii'));

DControlN=DControl-DDownStrong;
DControlN(DControlN<0)=0;
DDownStrongN=DDownStrong-DControl;
DDownStrongN(DDownStrongN<0)=0;

out.vol=cat(4,DControlN,DDownStrongN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'ControlDownStrongN.nii'));

ControlN=DControl-DDownStrong;
DControlN(DControlN<0)=0;
DDownStrongN=DDownStrong-DControl;
DDownStrongN(DDownStrongN<0)=0;

out.vol=cat(4,DControlN,DDownStrongN);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'ControlDownStrongN.nii'));