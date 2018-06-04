function f=SVDTailPoly(A,L)
for i=1:size(L,2)
%f(i)=A(1)+A(2)./L(i)+A(3)*exp(A(4)*L(i))+A(5).*(L(i).^-2);
f(i) = A(1)*exp(A(2)*L(i)) +A(3)+ A(4)*L(i);
end
end