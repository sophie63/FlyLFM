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
Data(Data<0)=0;

%Only keep the times when there is a clear flow in one direction compared
%to the other
 LeftminRight=Xk(:,1)-Xk(:,2);
 Left2=Xk(:,1);
 Left2(LeftminRight<60)=0;
 Right2=Xk(:,2);
 Right2((-LeftminRight)<60)=0;

%Left2=squeeze(Xk(:,1));
%Right2=squeeze(Xk(:,2));

%% Average data during walk, groom, rest and turn
Dleft=zeros(S(1),S(2),S(3),1);
Dright=zeros(S(1),S(2),S(3),1);

for i =1:S(4)
%for i =1000:4000    
    Dleft(:,:,:,1)=Dleft(:,:,:,1)+Data(:,:,:,i)*Left2(i);
    Dright(:,:,:,1)=Dright(:,:,:,1)+Data(:,:,:,i)*Right2(i);
end

Dleft=Dleft/(sum(Left2));
Dright=Dright/(sum(Right2));

 DleftrightN=Dleft-Dright;
 D1=DleftrightN;
 D1(D1<0)=0;
 D2=-DleftrightN;
 D2(D2<0)=0;
%D1=Dleft;
M=prctile(reshape(D1,S(1)*S(2)*S(3),1),99)
D1=D1/M;
M2=max(max(max(D1)));

%D2=Dright;
M=prctile(reshape(D2,S(1)*S(2)*S(3),1),99)
D2=D2/M;
M2b=max(max(max(D2)));
D1=D1/max(M2,M2b);
D2=D2/max(M2,M2b);
out.vol=cat(4,D1,D2);
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'LeftRight.nii'));


% It is unclear for me why the normalisation is needed/justified here. How
% come the intensity of left and right are not the same after normalisation
% with the amount of right or left?
% Could be problem with extraction left and right. Could be that the
% intensity of the activity is dependant on wether turn left or right, but
% not the speed. In that case thresholding and binarizingLeft and Right
% should give the same result.

% The advantage of this vs binarizing is that times when we are more
% confident about turning (high turning value) are weighted more

D1m=Montage3(D1);
Sm=size(D1m);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); 

fullFileName = fullfile(strcat(file(1:size(file,2)-4),'LeftRightM.PNG'));
imwrite(Dm, fullFileName);

