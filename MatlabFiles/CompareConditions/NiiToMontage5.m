[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

D1=Data(:,:,:,1);
D2=Data(:,:,:,2);

D1m=Montage5(D1);
D2m=Montage5(D2);
D3=zeros(size(D2));
D3m=Montage5(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
M=max(max(max(max(Dm4norm))));
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'OdorFlashM.PNG'));
imwrite(Dm/M, fullFileName);
