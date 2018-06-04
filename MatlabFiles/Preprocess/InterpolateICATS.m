clear all

% Open Time video
[FileName,PathName] = uigetfile('*.mat','Select the time for the video');
load(strcat(PathName,FileName))

% Open Time camera
[FileName,PathName] = uigetfile('*.mat','Select the time for the fluo');
load(strcat(PathName,FileName))

% Open Data
[FileName,PathName] = uigetfile('*.mat','Select the time series file');
file=strcat(PathName,FileName)
load(file)


%TimefluoICA=TimeFluoOn;
TimefluoICA=TimeFluoOn(2:(end-1));

T=TimeFluoOnVid';

for i=1:size(T,2)
    if T(i)>TimefluoICA(1,1)
        break
    end
end
size(T,2)
Initvid=i;

for i=1:size(T,2)
    if T(i)>TimefluoICA(end)
        break
    end
end
Endvid=i-1;

if Endvid==0
    Endvid=size(T,2)-1;
end

Initvid
Endvid

Tvid=T(Initvid:Endvid);




for i=1:size(TSo,2)
       TS2(:,i)=interp1(TimefluoICA',squeeze(TSo(:,i)),Tvid',[]);
end

save(strcat(file(1:size(file,2)-4),'int.mat'),'TS2')
