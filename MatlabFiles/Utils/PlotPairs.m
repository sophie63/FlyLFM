% For Fig. S24

fmask = fileread('/media/test15/ArclightCombo/FlipFlopComponents.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/test15/ArclightCombo/FlipflopCompnumber.txt');
Dlist = strsplit(fdata);


figure 
hold on
i=0
for idx=1:length(Mlist)

    file=Mlist{idx}
    file2=strcat(file(1:size(file,2)-6),'TS.mat');
    load(file2)
    ICnum=str2num(Dlist{idx})

    TS(:,1)=TSo(:,ICnum(1))/std(TSo(:,ICnum(1)));
    TS(:,2)=TSo(:,ICnum(2))/std(TSo(:,ICnum(2)));
    
    %ICmean=(TS(:,1)+TS(:,2))/2;
    plot(TS(:,1)- i*13)
    plot(TS(:,2)- i*13+5)
    
    %hold on
    %plot(TS(:,2)-ICmean+4)
    clear TSo
    clear TS
    i=i+1
end