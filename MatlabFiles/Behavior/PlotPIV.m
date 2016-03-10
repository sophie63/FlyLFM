%This routine opens the RGB avi resulting from optic flow analysis in FIJI

figure
clear

[FileName,PathName] = uigetfile('*.avi','Open colored PIV video','/media/');
file2=strcat(PathName,FileName)
uiopen(file2,1)

name=strcat('x',FileName)
[A,B]=strsplit(name,'.')
D=eval(A{1});

FlyLeft(1)=0;
FlyRight(1)=0;

for i=1:12116
FlyLeft(i+1)=mean(mean(D(i).cdata(:,:,2)));
FlyRight(i+1)=mean(mean(D(i).cdata(:,:,3)));
end

FlyLeft(FlyLeft>60)=60;
FlyRight(FlyRight>60)=60;

plot(smooth(FlyLeft),'.','color','g')
hold on 
plot(smooth(FlyRight),'.','color','b')
