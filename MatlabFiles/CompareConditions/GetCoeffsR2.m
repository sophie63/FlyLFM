
load('87to75.mat')
TS=TS';
for i=1:75
mdl = fitlm(Rkd(1:4179),TS(:,VarName3(i)));
R2(i)=mdl.Rsquared.Ordinary;
Coef(i)=mdl.Coefficients.Estimate(2);
end
Coef=Coef';
R2=R2';
