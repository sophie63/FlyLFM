
px1bin(1:8)=de2bi(67,8);
px2bin(1:8)=de2bi(31,8);
px3bin(1:8)=de2bi(87,8);
px4bin(1:8)=de2bi(52,8);

pxbin=[px4bin,px3bin,px2bin,px1bin];
sec=pxbin(26:32);
cyc=pxbin(13:25);
cycoffset=pxbin(1:12);

cycnum=bi2de(cyc)
secnum=bi2de(sec)
cycoffsetnum=bi2de(cycoffset)
T=double(secnum)+((double(cycnum)+double(cycoffsetnum)/3072.0)/8000.0);

j=0;

T2(1)=T(1);

for i=1:(size(T,2)-1)
Der=T(i+1)-T(i);
if Der<0
  j=j+1;
end
T2(i+1)=T(i+1)+j*128;
end

Time=T2-T2(1);
max(Time)

save(strcat(name(1:size(name,2)-4),'Time.mat'),'T')
