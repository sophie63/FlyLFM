clear all
%open mask
[FileName,PathName] = uigetfile('*.nii','Select the Mask files','MultiSelect','on');
files=strcat(PathName,FileName);
[FileName,PathName] = uigetfile('*.nii','Select the data files','MultiSelect','on');
filesData=strcat(PathName,FileName);
for j=1:size(files,2)
clear D Data S Data2 out R R2 S1 B D2 TS Av M2 Mask M
file=files{j};
M=MRIread(file);
Mask=M.vol;

%Mask=1-Mask;%use that line if mask from segmentation editor
%Md=double(M);

%open data
fileData=filesData{j};
D=MRIread(fileData);
Data=D.vol;

S=size(Data);
M2=Mask./max((max(max(max(Mask)))));
SM2=size(M2);
Data(isnan(Data))=0;

for i=1:S(3)
    for j=1:SM2(3)
        Av(i,j)=sum(sum(sum(M2(:,:,j).*squeeze(Data(:,:,i)))))/sum(sum(sum(M2(:,:,j))));
    end
end

save(strcat(file(1:size(file,2)-4),'AvROI.mat'),'Av')
csvwrite(strcat(file(1:size(file,2)-4),'TSROI.csv'),Av)
TS=Av';
plotTS
saveas(gcf,strcat(file(1:size(file,2)-4),'TimeSeries.png'))
close all
end
%load('87to75.mat')
%Avgood=Av(:,VarName3);
%Avgood=Avgood';