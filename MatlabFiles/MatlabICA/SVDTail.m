function f=SVDTail(sig,L,Ndim,Tdim)
q=Ndim*Tdim;
LM=sig*sig*(1+sqrt(1/q))^2
Lm=sig*sig*(1-sqrt(1/q))^2
for i=1:size(L,2)
f(i)=(q/(2*pi*sig*sig))*sqrt((LM-L(i))*(L(i)-Lm))/L(i);
end
end