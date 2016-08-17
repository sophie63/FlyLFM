

function f=SVDTail(sig,L)
q=Ndim*Tdim;
LM=sig*sig*(1+sqrt(1/q))^2;
Lm=sig*sig*(1-sqrt(1/q))^2;
f=(q/(2*pi*sig*sig))*sqrt((LM-L)*(L-Lm))/L;
end

sig0 = 1;
Fsumsquares = @(x)sum((SVDTail(sig,L) - y).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[xunc,ressquared,eflag,outputu] = ...
    fminunc(Fsumsquares,sig0,opts)
