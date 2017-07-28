clear all

[FileName,PathName] = uigetfile('*.avi','Select the AVI file');
name=strcat(PathName,FileName)

mov=VideoReader(name);
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

save(strcat(name(1:size(name,2)-4),'Time.mat'),'Time')
