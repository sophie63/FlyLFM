clear all

% Open Time video
[FileName1,PathName1] = uigetfile('*.mat','Select the time for the video');


% Open Time camera
[FileName2,PathName2] = uigetfile('*.mat','Select the time for the fluo');


% Open Data
[FileName3,PathName3] = uigetfile('*.nii','Select the Nifti file');


file=strcat(PathName3,FileName3)
B=MRIread(file);
D=double(B.vol);
clear B
load(strcat(PathName1,FileName1))
load(strcat(PathName2,FileName2))

TimefluoICA=TimeFluoOn(2:(end-1))';

%T=TimeFluoOnVid;
T=Time/1000;

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
S=size(D);
for i=1:S(1)
    for k=1:S(3)
        parfor l=1:S(2)
            Data2(i,l,k,:)=interp1(TimefluoICA',squeeze(D(i,l,k,:)),Tvid',[]);
        end
    end
    i
end

out.vol=Data2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'int.nii'));

Initvid
Endvid
