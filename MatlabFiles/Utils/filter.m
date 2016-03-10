function [D2]=filter(D)

S=size(D)

for i=1:S(1)
    for j=1:S(2)
        for k=1:S(3)
            D2(i,j,k,:)=smooth(D(i,j,k,:),2);
        end
    end
end
