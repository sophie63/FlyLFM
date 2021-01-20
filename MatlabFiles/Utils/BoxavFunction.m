function Ds=BoxavFunction(Data)

S=size(Data);
d=3;
Ds=zeros(floor(S(1)/d), floor(S(2)/d),S(3)/d,S(4));

for j=1:S(4)
        X=squeeze(Data(1:floor(S(1)/d)*d,1:floor(S(2)/d)*d,1:floor(S(2)/d)*d,j));
        Xs=reshape(sum(sum(reshape(X, d, floor(S(1)/d), d, floor(S(2)/d)), d, floor(S(3)/d))), 3), floor(S(1)/d), floor(S(2)/d,floor(S(3)/d));
        Ds(:,:,:,j)=Xs;
end
    
end


