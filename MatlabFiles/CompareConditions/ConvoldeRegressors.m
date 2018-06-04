I=[192:4772];

%X=zeros(size(Left,2),6);
X(:,1)=Left(1,I);
X(:,2)=Right(1,I);
% X(:,1)=Left;
% X(:,2)=Right;
%X(:,3)=Walk;
%X(:,5)=Groom;
%X(:,6)=Rest;
%X(:,3)=Straight;
%X(:,4)=WalkHand;
%X(:,5)=Groom;
%Fr=60;
Fr=50;
S=size(X)
for i=1:S(2)
    Xk2(:,i)=conv(GCamp6Fkernel(1:100/Fr:202),X(:,i));
end
Xk=Xk2(1:S(1),:);