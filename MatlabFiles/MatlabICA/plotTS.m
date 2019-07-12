
f=figure

clear C

 C(:,1)=[1,0,0];
 C(:,2)=[0,1,0];
 C(:,3)=[0,0,1];
 C(:,4)=[0.8,0.8,0];
 C(:,5)=[0,1,1];
 C(:,6)=[1,0,1];

D=TS';
 
Imax=size(D,1);
h=10;


for i=1:size(D,2)
i-ceil((i-6)/6)*6
plot((D(:,i)-mean(D(:,i)))/sqrt(var(D(:,i)))+h*i,'color',C(:,i-ceil((i-6)/6)*6));
x1 = 0.005*Imax+2;
y1 = h*i;
str1 = int2str(i);
%text(x1,y1,str1)
hold on
end

%savefig('RegionTS.fig')
%ylim([-h h*(i+10)])
%  plot(Tvid(503:3557)-Tvid(503),13-Walk862*5+h*i,'r+')
%  plot(Tvid(503:3557)-Tvid(503),13-Idle862*5+h*i,'b+')
%  plot(Tstim,Stim-2,'m+')
%  plot(Tstim,Stim + h*i + 15,'m+')
