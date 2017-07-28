function [Straight,Left, Right] = OpticFlow2Direction(B)



S=size(B)

Right(1)=(4/3)*mean(mean(B(1).cdata(:,:,3)));
Left(1)=(4/3)*mean(mean(B(1).cdata(:,:,2)));
Straight(1)=-(1/4)*Right(1)-(1/4)*Left(1)+mean(mean(B(1).cdata(:,:,1)));

for i=1:(S(2)-1)
Right(i+1)=0.5*(4/3)*((mean(mean(B(i).cdata(:,:,3)))+mean(mean(B(i+1).cdata(:,:,3)))));
Left(i+1)=0.5*(4/3)*(mean(mean(B(i).cdata(:,:,2)))+mean(mean(B(i).cdata(:,:,2))));
Straight(i+1)=-(1/4)*Right(i)-(1/4)*Left(i)+0.5*(mean(mean(B(i).cdata(:,:,1)))+mean(mean(B(i).cdata(:,:,1))));
end

Right(S(2)+1)=(4/3)*mean(mean(B(S(2)).cdata(:,:,3)));
Left(S(2)+1)=(4/3)*mean(mean(B(S(2)).cdata(:,:,2)));
Straight(S(2)+1)=-(1/4)*Right(S(2)+1)-(1/4)*Left(S(2)+1)+mean(mean(B(S(2)).cdata(:,:,1)));


