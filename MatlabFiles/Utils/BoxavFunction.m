function Ds=BoxavFunction(Data)

S=size(Data);

Ds=zeros(floor(S(1)/5), floor(S(2)/5),S(3)/5,S(4));

for j=1:S(4)
        X=squeeze(Data(1:floor(S(1)/5)*5,1:floor(S(2)/5)*5,1:floor(S(2)/5)*5,j));
        Xs=reshape(sum(sum(reshape(X, 5, floor(S(1)/5), 5, floor(S(2)/5)), 5, floor(S(3)/5))), 3), floor(S(1)/5), floor(S(2)/5,floor(S(3)/5));
        Ds(:,:,:,j)=Xs;
end
    
end


