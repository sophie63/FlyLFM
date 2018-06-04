figure

params.Fs=50;
params.pad=2;

[mspectrumCXGCaMP6(:,1),f]=mtspectrumc(diff(TSGCaMP6(:,8)),params);
[mspectrumCXGCaMP6(:,2),f]=mtspectrumc(diff(TSGCaMP6_2(:,8)),params);
[mspectrumCXGCaMP6(:,3),f]=mtspectrumc(diff(TSGCaMP6_3(:,8)),params);
[mspectrumCXGCaMP6(:,4),f]=mtspectrumc(diff(TSGCaMP6_4(:,8)),params);
params.Fs=200;
SpectrumGCaMP6=mean(mspectrumCXGCaMP6,2);
plot_vector(SpectrumGCaMP6,f,'l')
hold on
[mspectrumCXArcLight(:,1),f]=mtspectrumc(diff(TSArclight_2(:,8)),params);
[mspectrumCXArcLight(:,2),f]=mtspectrumc(diff(TSArclight_3(:,8)),params);
[mspectrumCXArcLight(:,3),f]=mtspectrumc(diff(TSArclight_4(:,8)),params);
[mspectrumCXArcLight(:,4),f]=mtspectrumc(diff(TSArclight_5(:,8)),params);
[mspectrumCXArcLight(:,5),f]=mtspectrumc(diff(TSArclight_6(:,8)),params);
SpectrumArclight=mean(mspectrumCXArcLight,2);

figure
plot_vector(SpectrumArclight,f,'l')

params.Fs=200;
[mspectrumCXcontrol(:,1),f]=mtspectrumc(diff(TScontrol(:,8)),params);
figure
plot_vector(mspectrumCXcontrol,f,'l')
params.Fs=100;
[mspectrumCXcontrol2,f]=mtspectrumc(diff(TScontrol_2(500:end,8)),params);
figure
plot_vector(mspectrumCXcontrol2,f,'l')