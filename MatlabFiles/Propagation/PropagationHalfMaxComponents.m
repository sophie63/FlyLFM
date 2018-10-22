clear all

prompt = 'What is the time of onset?';
To = input(prompt)
% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);

[FileName,PathName] = uigetfile('*.mat','Select the region time series');
file2=strcat(PathName,FileName)
load(file2)
%load(strcat(file(1:size(file,2)-6),'TS.mat'))


D1=zeros(S1(1),S1(2),S1(3),1);
D2=zeros(S1(1),S1(2),S1(3),1);

for i=1:S1(4)
    Before=mean(TSo(To-50:To,i));
    StdBefore=std(TSo(To-50:To,i));
    MaxResponse=max(TSo(To:To+60,i));
    HalfMax=Before+(MaxResponse-Before)/2;
    if (MaxResponse-Before)>2*StdBefore
        [c index]=min(abs(TSo(To:(To+40),i)-HalfMax));
        Delay1(i)=index;
        D1=D1+D(:,:,:,i)*(40-Delay1(i))/40;
        D2=D2+D(:,:,:,i)*Delay1(i)/40;
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
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'DelaysHalfMax.PNG'));
imwrite(Dm/M, fullFileName);



D1=zeros(S1(1),S1(2),S1(3),1);
D2=zeros(S1(1),S1(2),S1(3),1);

for i=1:S1(4)
    Before=mean(TSo(To-50:To,i));
    StdBefore=std(TSo(To-50:To,i));
    MaxResponse=max(TSo(To:To+50,i));
    HalfMax=Before+(MaxResponse-Before)/2;
    if (MaxResponse-Before)>2*StdBefore
        [c index]=min(abs(TSo(To:(To+70),i)-MaxResponse));
        Delay2(i)=index;
        D1=D1+D(:,:,:,i)*(70-Delay2(i))/70;
        D2=D2+D(:,:,:,i)*Delay2(i)/70;
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
fullFileName = fullfile(strcat(file(1:size(file,2)-4),'DelaysMax.PNG'));
imwrite(Dm/M, fullFileName);








