function A=bit2num(B)

A=0;
S=size(B);

for i=1:S(2)
    A=A+double(B(S(2)-i+1)*(2^(i-1)));
end
