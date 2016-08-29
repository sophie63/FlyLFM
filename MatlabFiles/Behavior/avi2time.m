function Time = avi2time (name)

mov=VideoReader(name);
B = read(mov,[1 Inf]);
S=size(B)
for i=1:S(4)
A=B(:,:,:,i);
D(:,i)=A(1,1:4,1); 
Px1(i)=int64(D(1,i));
Px2(i)=int64(D(2,i));
Px3(i)=int64(D(3,i));
px2bin(i,:)=byteestobit(Px2(i));
px1bin(i,:)=byteestobit(Px1(i));
px3bin(i,:)=byteestobit(Px3(i));
end

pxbin=[px1bin,px2bin,px3bin];
sec=pxbin(:,1:7);
cyc=pxbin(:,8:20);

for i=1:S(4)
cycnum(i)=bit2num(cyc(i,:));
end

for i=1:S(4)
secnum(i)=bit2num(sec(i,:));
end

Time=secnum+cycnum*0.000125;

j=0;

T2(1)=T(1);

for i=1:(size(T,2)-1)
Der=T(i+1)-T(i);
if Der<0
   j=j+1
   i
end
T2(i+1)=T(i+1)+j*128;
end

Time=T2-T2(1);
