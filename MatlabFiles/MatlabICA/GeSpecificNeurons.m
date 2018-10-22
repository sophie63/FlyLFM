load('87to75.mat')

Fs=200;

Sm2=Sm(:,VarName3);
%Sm2=Av(VarName3,:);
%Sm2=Sm2';
S=size(Sm2,1);
j=0;
figure 
hold on

SMbin=zeros(size(Sm2));
SMbin(Sm2>0)=1;

for i=1:S
    SortedSm=sort(Sm2(i,:));
  if (Sm2(i,22)>0 && Sm2(i,6)>0 && Sm2(i,34)>0) && (Sm2(i,22)==max(Sm2(i,:)) || Sm2(i,6)==max(Sm2(i,:)) || Sm2(i,34)==max(Sm2(i,:)))
 %    if (Sm2(i,22)==max(Sm2(i,:)) || Sm2(i,6)==max(Sm2(i,:)) || Sm2(i,34)==max(Sm2(i,:))) && sum(SMbin(i,:))<15 && (Sm2(i,21)>5*(SortedSm(4)))
     
       'AL_PN right' 
       i
       plot(TSo(6500:end,i)+0.5*j)
       j=j+1;      
       MedFreq((TSo(6500:end,i)),Fs)
       ALPN{2}{end+1}=i;
       ALPN{3}{end+1}=MedFreq((TSo(6500:end,i)),Fs);
       X=[ones(size(Xk(:,2))) Xk(:,2)];
       [b,bint,r,rint,stats]=regress(TSo(:,i),X);
       ALPN{4}{end+1}=sqrt(stats(1));
   end
end

j=0;
figure 
hold on
for i=1:S
    SortedSm=sort(Sm2(i,:));
   if (Sm2(i,45)>0 && Sm2(i,59)>0 && Sm2(i,70)>0) && (Sm2(i,45)==max(Sm2(i,:)) || Sm2(i,59)==max(Sm2(i,:)) || Sm2(i,70)==max(Sm2(i,:)))
   %  if (Sm2(i,45)==max(Sm2(i,:)) || Sm2(i,59)==max(Sm2(i,:)) || Sm2(i,70)==max(Sm2(i,:))) && sum(SMbin(i,:))<15 && (Sm2(i,21)>5*(SortedSm(4)))

    'AL_PN left'
       i
       plot(TSo(6500:end,i)+0.5*j)
       j=j+1;
       MedFreq((TSo(6500:end,i)),Fs)
       ALPN{2}{end+1}=i;
       ALPN{3}{end+1}=MedFreq((TSo(6500:end,i)),Fs);
       X=[ones(size(Xk(:,2))) Xk(:,2)];
       [b,bint,r,rint,stats]=regress(TSo(:,i),X);
       ALPN{4}{end+1}=sqrt(stats(1));
   end
end

j=0;
figure 
hold on


for i=1:S
   SortedSm=sort(Sm2(i,:));
   if (Sm2(i,21)==max(Sm2(i,:))) && (sum(SMbin(i,:))<30) && (Sm2(i,21)>5*(SortedSm(4)))
        'EB'
        i
        j=j+1;
        MedFreq((TSo(:,i)),Fs)
   end
end

for i=1:S
   SortedSm=sort(Sm2(i,:));
   if (Sm2(i,5)==max(Sm2(i,:))) && (sum(SMbin(i,:))<30) && (Sm2(i,5)>5*(SortedSm(4)))
        'PB'
        i
        j=j+1;
        MedFreq((TSo(:,i)),Fs)
   end
end

