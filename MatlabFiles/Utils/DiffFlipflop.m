% For Fig. S24

fmask = fileread('/media/test4/ArclightCombo/FlipFlopComponents.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/test4/ArclightCombo/FlipflopCompnumber.txt');
Dlist = strsplit(fdata);


figure 
hold on
i=0
for idx=1:length(Mlist)
    figure
    file=Mlist{idx}
    file2=strcat(file(1:size(file,2)-6),'TS.mat');
    load(file2)
    ICnum=str2num(Dlist{idx})

    TS(:,1)=TSo(:,ICnum(1))/std(TSo(:,ICnum(1)));
    TS(:,2)=TSo(:,ICnum(2))/std(TSo(:,ICnum(2)));
    
    ICmean=(TS(:,1)+TS(:,2))/2;
    for i=1:2
        plot((1:size(TS,1))/200,TS(:,i)+(i-1)*4)
  
        hold on
    end
    
    set(gca,'ytick',[])
    xlabel('Time (s)')
    saveas(gcf,strcat('Fig',num2str(idx),'.svg'))
    %hold on
    %plot(TS(:,2)-ICmean+4)
    clear TSo
    clear TS
    i=i+1
end

