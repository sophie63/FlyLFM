function Ds=BoxavFunction3D(Data)

S=size(Data);

Ds=zeros(floor(S(1)/5), floor(S(2)/5),floor(S(3)/5));

X=squeeze(Data(1:floor(S(1)/5)*5,1:floor(S(2)/5)*5,1:floor(S(3)/5)*5));
size(X)
Xs=reshape(sum(sum(sum(reshape(X, 5, floor(S(1)/5), 5, floor(S(2)/5), 5, floor(S(3)/5)),1),3),5), floor(S(1)/5), floor(S(2)/5),floor(S(3)/5));
Ds(:,:,:)=Xs;
   
end


