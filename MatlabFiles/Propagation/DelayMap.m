
% Choose the reconstructed 4D nifti file where all values are positive
%(result of the addmin program)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

Tflash=13698;
Todor=3888;
FR=50;

Sts=size(TSo);
DelayFlash=zeros(1,S(2));

TSfb=mean(TSo((Tflash-FR):Tflash,:),1);
TSfa=mean(TSo(Tflash:(Tflash+FR),:),1);
TSfstd=std(TSo((Tflash-FR):Tflash,:),1);

D1=zeros(S(1),S(2),S(3),1);
D2=zeros(S(1),S(2),S(3),1);
for i=1:Sts(2)
    if (TSfa(i)-TSfb(i))>6*TSfstd(i)
        M=max(TSo(Tflash:(Tflash+FR),i)-TSfb(i))/2;
        [c index]=min(abs(TSo(Tflash:(Tflash+FR),i)-TSfb(i)-M));
        DelayFlash(i)=index;
        D1=D1+Data(:,:,:,i)*(FR-DelayFlash(i))/FR;
        D2=D2+Data(:,:,:,i)*DelayFlash(i)/FR;
    end
end



D1m=Montage3(D1);
D2m=Montage3(D2);
D3=zeros(size(D2));
D3m=Montage3(D3);
Dm=cat(3,D1m,D2m,D3m);
Dm4norm=Dm;
Dm4norm(Dm==1)=0;
Sn=size(Dm4norm); 
M=prctile(reshape(Dm4norm,Sn(1)*Sn(2)*Sn(3),1),99.9);
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'DelaysAfterFlashM.PNG'));
imwrite(Dm/M, fullFileName);





