params.Fs=50;
params.pad=2;

[mspectrum,f]=mtspectrumc(diff(TS(:,4)),params);
figure
plot_vector(mspectrum,f,'l')

movingwin=[5 1];
[S,t,f]=mtspecgramc(diff(TS(:,4)),movingwin,params);
figure
plot_matrix(S,t,f)