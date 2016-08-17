sig0 = 1;
Fsumsquares = @(sig)sum((SVDTail(sig,500:1000,S1(1)*S1(2)*S1(3),S1(4)) - diag(ss(500:1000,500:1000))').^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[sigunc,ressquared,eflag,outputu] = fminunc(Fsumsquares,sig0,opts);
