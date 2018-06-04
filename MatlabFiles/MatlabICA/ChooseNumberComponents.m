Dvalue=max(Spectrum)-min(Spectrum);
Dvalue=Dvalue/1000;
plot(smooth(diff(smooth(Spectrum,20))))