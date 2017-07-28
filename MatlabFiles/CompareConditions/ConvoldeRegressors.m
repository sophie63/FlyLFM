I=[(270+7):(11174+270)];

%X(:,1)=Left(I);
%X(:,2)=Right(I);
%X(:,3)=Straight(I);
%X(:,4)=WalkHand;
%X(:,5)=Groom;
Fr=50;

S=size(X)
for i=1:S(2)
    Xk(:,i)=conv(GCamp6Fkernel(1:100/Fr:202),X(:,i));
end
Xk=Xk(1:S(1),:);