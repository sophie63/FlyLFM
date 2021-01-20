function Ds=BoxavFunction3D(Data)

S=size(Data);
d=3;
Ds=zeros(floor(S(1)/d), floor(S(2)/d),floor(S(3)/5));

X=squeeze(Data(1:floor(S(1)/d)*d,1:floor(S(2)/d)*d,1:floor(S(3)/5)*5));
size(X)
Xs=reshape(mean(mean(mean(reshape(X, d, floor(S(1)/d), d, floor(S(2)/d), 5, floor(S(3)/5)),1),3),5), floor(S(1)/d), floor(S(2)/d),floor(S(3)/5));
Ds(:,:,:)=Xs;
   
end


