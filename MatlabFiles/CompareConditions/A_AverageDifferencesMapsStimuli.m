clear all


prompt = 'What is the time of the first odor?';
Todor = input(prompt)

prompt = 'What is the time of the first flash?';
Tflash = input(prompt)

prompt = 'What is the frame rate?';
Fr = input(prompt)


% Choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);


%% Average data during walk, groom, rest and turn
Dodor=zeros(S(1),S(2),S(3),1);
Dflash=zeros(S(1),S(2),S(3),1);

Dodor=mean(Data(:,:,:,Todor:Todor+Fr),4)-mean(Data(:,:,:,Todor-Fr:Todor),4);
Dflash=mean(Data(:,:,:,Tflash:Tflash+Fr),4)-mean(Data(:,:,:,Tflash-Fr:Tflash),4);

%Ddiff=Dodor-Dflash;
%D1=DleftrightN;
%D1(DleftrightN<0)=0;
%D2=-DleftrightN;
%D2(DleftrightN>0)=0;
D1=Dodor/max(max(max(Dodor)));
D2=Dflash/max(max(max(Dflash)));
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'OdorFlash.nii'));

D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'OdorFlashM.PNG'));
imwrite(Dm/M, fullFileName);
