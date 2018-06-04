
Fsumsquares = @(A) sum((SVDTailPoly(A,2:500)-Spectrum(2:500)').^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
A0=[6, -0.02, 6, -0.005];
[A,ressquared,eflag,outputu] = fminunc(Fsumsquares,A0,opts);