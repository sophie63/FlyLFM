function Dmont = Montage3(D)

% First project over 3 layers
S1=size(D);
Inter=floor(S1(3)/5);
for l=1:5
   D3(:,:,l)=mean(D(:,:,(l-1)*Inter+1:l*Inter),3);
end

% Combine
S=size(D3);


D3b=ones(S(1)+1,S(2),S(3));
D3b(1:S(1),:,:)=D3;

Dmont=cat(1,D3b(:,:,1),D3b(:,:,2),D3b(:,:,3),D3b(:,:,4),D3b(:,:,5));
Dmont=Dmont(1:S(1)*5+2,:);
