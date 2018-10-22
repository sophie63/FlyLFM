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
    
    ICmean=(TS(:,1)+TS(:,2))/2;
    plot(TS(:,1)-ICmean + i*3)
    %hold on
    %plot(TS(:,2)-ICmean+4)
    clear TSo
    clear TS
    i=i+1
end

load('/media/test15/ControlCombo/100353/Fliflop/100353ss2oncregcdFF20sTSROI.mat')
TS(:,1)=TS(:,1)/std(TS(:,1));
TS(:,2)=TS(:,2)/std(TS(:,2));
ICmean=(TS(:,1)+TS(:,2))/2;
plot(TS(:,1)-ICmean-6)
load('/media/test15/ControlCombo/100682series/FlipFlop/100682ss1ondFF20sTSROI.mat')
TS(:,1)=TS(:,1)/std(TS(:,1));
TS(:,2)=TS(:,2)/std(TS(:,2));
ICmean=(TS(:,1)+TS(:,2))/2;
plot(TS(:,1)-ICmean-9)