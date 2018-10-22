params.Fs=50;
params.pad=2;

[mspectrum,f]=mtspectrumc(diff(TSo(:,633)),params);
figure
plot_vector(mspectrum,f,'l')

movingwin=[5 1];
[S,t,f]=mtspecgramc(diff(TSo(:,633)),movingwin,params);
figure
plot_matrix(S,t,f)

s=diff(TSo(:,633));

Fs = 50;                                   % Sampling Frequency
Fn = Fs/2;                                              % Nyquist Frequency
L = length(TSo(:,633));
FTs = fft(s)/L;
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                     % Frequency Vector
Iv = 1:length(Fv);  
CumAmp = cumtrapz(Fv, abs(FTs(Iv)));                    % Integrate FFT Amplitude
MedFreq = interp1(CumAmp, Fv, CumAmp(end)/2);  
figure(1)
plot(Fv, abs(FTs(Iv))*2, '-b')                          % Plot FFT
hold on
plot(Fv, CumAmp, '-g')                                  % Plot Cumulative Amplitude Integral
plot([MedFreq MedFreq], ylim, '-r', 'LineWidth',1)      % Plot Median Frequency
hold off
grid
