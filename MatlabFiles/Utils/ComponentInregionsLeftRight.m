fmask = fileread('/media/test5/ComponentsList.txt');
Mlist = strsplit(fmask);

load('/media/test5/FreeBehaviorPanNeuronalGCaMP6/Components_cbyr.mat')
Components_cbyr(Components_cbyr==0)=NaN;

i=1
for idx=1:length(Mlist)
   
    file=Mlist{idx}    
    file2=strcat(file(1:size(file,2)-4),'ICRegions.mat')
    load(file2)
    S=size(Av);
    for j=1:4
        if isnan(Components_cbyr(i,j))
            Avcbyr(:,j,i)=NaN;
        else    
        Avcbyr(:,j,i)=Av(:,int8(Components_cbyr(i,j)-1));
        end
    end
    
    i=1+i
    clear Av
end

Avcb=zeros(87,2,7);
Avcb(:,1,:)=nanmean(Avcbyr(:,[1,3],:),2);
Avcb(:,2,:)=nanmean(Avcbyr(:,[2,4],:),2);
load('87to75.mat')
Avcbgood=Avcb(VarName3,:,:);



load('Regions75to41.mat')
for i=1:41
if Region75to41(i,1)==0
Avcb41(i,:,:)=Avcbgood(i,:,:);
else
Avcb41(i,:,:)=mean(Avcbgood(Region75to41(i),:,:),1);
end
end

AvcbPan=nanmean(Avcb41(:,:,[1,2,3,5,6,7]),3);
