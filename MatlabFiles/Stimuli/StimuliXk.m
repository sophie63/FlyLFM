% Tstim=np.arange(0,np.max(TimeFluoICA),0.001)
% Tstim2=Tstim-np.max(Time_fluo862)
% Flashes=np.zeros(Tstim.shape[0])
% Odor=np.zeros(Tstim.shape[0])
% for i in range(Tstim.shape[0]):
%     if (Tstim[i]>10 and Tstim[i]<12) or (Tstim[i]>14 and Tstim[i]<16) or (Tstim[i]>18 and Tstim[i]<20) or (Tstim[i]>22 and Tstim[i]<24)or (Tstim[i]>26 and Tstim[i]<28)or (Tstim[i]>30and Tstim[i]<32):
%         Flashes[i]=1
%     if (Tstim2[i]>12 and Tstim2[i]<14) or (Tstim2[i]>19 and Tstim2[i]<21) or (Tstim2[i]>26 and Tstim2[i]<28):
%         Odor[i]=1
% Flashes[Flashes==0]=np.nan
% Odor[Odor==0]=np.nan 

FlashesOn=zeros(size(TimeFluoICA));
OdorOn=zeros(size(TimeFluoICA));

[~,I]=min(abs(TimeFluoICA-10));
FlashesOn(I)=1;
[~,I]=min(abs(TimeFluoICA-14));
FlashesOn(I)=1;
[~,I]=min(abs(TimeFluoICA-18));
FlashesOn(I)=1;
[~,I]=min(abs(TimeFluoICA-22));
FlashesOn(I)=1;
[~,I]=min(abs(TimeFluoICA-26));
FlashesOn(I)=1;
[~,I]=min(abs(TimeFluoICA-30));
FlashesOn(I)=1;

Toodor=35.595852098795213;
[~,I]=min(abs(TimeFluoICA-Toodor-12-0.35));
OdorOn(I)=1;
[~,I]=min(abs(TimeFluoICA-Toodor-19-0.35));
OdorOn(I)=1;
[~,I]=min(abs(TimeFluoICA-Toodor-26-0.35));
OdorOn(I)=1;


