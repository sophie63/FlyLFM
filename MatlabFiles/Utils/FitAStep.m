alternativeEquation = 'a+c*heaviside(x-p)*(x-p)';
disp(alternativeEquation)

st_ = [0 0 0];
ft2 = fittype(alternativeEquation, ...
    'dependent',{'y'},'independent',{'x'}, ...
    'coefficients',{'a','c','p'});
x=160:280;
x=x';
a0=mean(TS(2,160:280));

fo2 = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[-1 -1 200],...
    'Upper',[1 1 300], ...
    'Startpoint',[a0 0.0014 225]);

sample=TS(2,160:280)';
[cf2, gof2] = fit(x,sample,ft2,fo2);

clf
plot(x, sample);
hold on
plot(x,cf2(x),'g','linewidth',2);
set(gca,'Ylim',[-.5 .5])
legend({'data','fitted curve'})

cf2.p



for i=i:12
    sample=TS(i,160:550)';
    av=mean(sample(1:60));
    M=max(sample);
    [c index]=min(abs(sample-av-(M-av)/2));
    DelayFlash(i)=index;
end
    