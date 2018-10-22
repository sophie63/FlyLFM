function MFreq=MedFreq(s,Fs)

%s=diff(TSo(:,319));

%Fs = 50;                                   % Sampling Frequency
Fn = Fs/2;                                              % Nyquist Frequency
L = length(s);
FTs = fft(s)/L;
S=size(FTs);
%FTs = sgolayfilt(FTs,3,149);



for k=1:2
for i=1:(size(FTs,1)-10)
    if abs(FTs(i))>(mean(abs(FTs((i+1):(i+10))))+3*std(abs(FTs(i:(i+10)))))
        FTs(i)=FTs(i)/20;
    end
end
end



Fvfull = linspace(0, 1, fix(L/2)+1)*Fn;                     % Frequency Vector
% Use only the spectrum up to 20Hz
[C,Ind20Hz]=min(abs(Fvfull-20));
Fv = Fvfull(1:Ind20Hz);
Iv = 1:length(Fv);  

CumAmp = cumtrapz(Fv, abs(FTs(Iv)));   % Integrate FFT Amplitude


MFreq = interp1(CumAmp', Fv, CumAmp(end)/2);  

figure
plot(Fv, smooth(abs(FTs(Iv))*2), '-b')                          % Plot FFT
hold on
plot(Fv, CumAmp, '-g')                                  % Plot Cumulative Amplitude Integral
plot([MFreq MFreq], ylim, '-r', 'LineWidth',1)      % Plot Median Frequency
hold off
grid
