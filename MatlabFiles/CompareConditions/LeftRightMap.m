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

% LeftminRight=Xk(:,1)-Xk(:,2);
% Left2=LeftminRight;
% Left2(LeftminRight<0)=0;
% Right2=-LeftminRight;
% Right2(LeftminRight>0)=0;

Left2=squeeze(Xk(:,1));
Right2=squeeze(Xk(:,2));

%% Average data during walk, groom, rest and turn
Dleft=zeros(S(1),S(2),S(3),1);
Dright=zeros(S(1),S(2),S(3),1);

for i =1:S(4)
    Dleft(:,:,:,1)=Dleft(:,:,:,1)+Data(:,:,:,i)*Left2(i);
    Dright(:,:,:,1)=Dright(:,:,:,1)+Data(:,:,:,i)*Right2(i);
end

Dleft=Dleft/(sum(Left2));
Dright=Dright/(sum(Right2));

% DleftrightN=Dleft-Dright;
% D1=DleftrightN;
% D1(DleftrightN<0)=0;
% D2=-DleftrightN;
% D2(DleftrightN>0)=0;
D1=Dleft;
M=prctile(reshape(D1,S(1)*S(2)*S(3),1),99.9);
D1=D1/M;
D2=Dright;
M=prctile(reshape(D2,S(1)*S(2)*S(3),1),99.9);
D2=D2/M;
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'LeftMinRight.nii'));

% It is unclear for me why the normalisation is needed/justified here. How
% come the intensity of left and right are not the same after normalisation
% with the amount of right or left?
% Could be problem with extraction left and right. Could be that the
% intensity of the activity is dependant on wether turn left or right, but
% not the speed. In that case thresholding and binarizingLeft and Right
% should give the same result.
D1m=Montage3(D1);
Sm=size(D1m);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); 

fullFileName = fullfile(strcat(file(1:size(file,2)-4),'LeftMinRightM.PNG'));
imwrite(Dm, fullFileName);
