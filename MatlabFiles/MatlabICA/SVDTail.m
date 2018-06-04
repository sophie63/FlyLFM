function f=SVDTail(sig,L,Ndim,Tdim)
q=Tdim;
LM=sig*sqrt(1+Ndim/q+2*sqrt(Ndim/q));
Lm=sig*sqrt(1+Ndim/q-2*sqrt(Ndim/q));
for i=1:size(L,2)
f(i)=(sqrt(q)/(2*pi*sig))*sqrt((LM*LM-L(i)*L(i))*(L(i)*L(i)-Lm*Lm));
end
end