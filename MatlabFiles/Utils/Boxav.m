for i=1:35
    for j=1:6262
        X=squeeze(D(:,:,i,j));
        Xs=reshape(sum(sum(reshape(X, 2, 42, 2, 87), 1), 3), 42, 87);
        Ds(:,:,i,j)=Xs;
    end
    i
end
