function Time = avi2time (B)

S=size(B)

for i=1:S(2)
A=B(i).cdata;
D(:,i)=A(1,1:4,1); 
Px1(i)=int64(D(1,i));
Px2(i)=int64(D(2,i));Px3(i)=int64(D(3,i));
px2bin(i,:)=byteestobit(Px2(i));
px1bin(i,:)=byteestobit(Px1(i));
px3bin(i,:)=byteestobit(Px3(i));
end

pxbin=[px1bin,px2bin,px3bin];
sec=pxbin(:,1:7);
cyc=pxbin(:,8:20);

for i=1:S(2)
cycnum(i)=bit2num(cyc(i,:));
end

for i=1:S(2)
secnum(i)=bit2num(sec(i,:));
end

Time=secnum+cycnum*0.000125;