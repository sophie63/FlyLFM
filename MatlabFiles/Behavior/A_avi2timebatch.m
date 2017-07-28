clear all

%choose the reconstructed 4D nifti file
[FileName,PathName] = uigetfile('*.avi','Select the AVI files','MultiSelect','on');
files=strcat(PathName,FileName);

for j=1:size(files,2)

file=files{j};
mov=VideoReader(file);
B = read(mov,[1 Inf]);
S=size(B)

for i=1:S(4)
A=B(:,:,:,i);
px1bin(i,1:8)=de2bi(squeeze(A(1,1,1)),8);
px2bin(i,1:8)=de2bi(squeeze(A(1,2,1)),8);
px3bin(i,1:8)=de2bi(squeeze(A(1,3,1)),8);
px4bin(i,1:8)=de2bi(squeeze(A(1,4,1)),8);
end

pxbin=[px4bin,px3bin,px2bin,px1bin];
sec=double(pxbin(:,26:32));
cyc=double(pxbin(:,13:25));
cycoffset=double(pxbin(:,1:12));

for i=1:S(4)
cycnum(i)=bi2de(cyc(i,:))*1000;
secnum(i)=bi2de(sec(i,:))*1000;
cycoffsetnum(i)=bi2de(cycoffset(i,:))*1000;
end

T=double(secnum)+((double(cycnum)+double(cycoffsetnum)/3072.0)/8000.0);

j=0;

T2(1)=T(1);

for i=1:(size(T,2)-1)
Der=T(i+1)-T(i);
if Der<0
  j=j+1;
end
T2(i+1)=T(i+1)+j*128000;
end

Time=T2-T2(1);
max(Time)

save(strcat(file(1:size(file,2)-4),'Time.mat'),'Time')

clear B A px1bin px2bin px3bin px4bin Time T2 pxbin sec cyc cycoffset S cycnum secnum cycoffsetnum T

end

