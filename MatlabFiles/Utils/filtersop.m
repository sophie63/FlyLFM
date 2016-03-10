function [D2]=filtersop(D)

S=size(D)

for i=1:S(1)
    for j=1:S(2)
        for k=1:S(3)
            D2(i,j,k,:)=medfilt1(D(i,j,k,:),3);
        end
    end
end
