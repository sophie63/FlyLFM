function [Trimed]=Trim2D(Trim,n)

% This function removes n pixels around the binary 2D shape

S=size(Trim);
Trim(Trim<0.5)=0;

Trima=Trim;
Trimint=Trim;
Trimb=Trim;
Trimintb=Trim;
for k=1:n
for i =1:S(1)
    for j= 1:(S(2)-1)
        if Trima(i,j)==0 && Trima(i,j+1)>0
            Trimint(i,j+1)=0;
        else if Trima(i,j)>0 && Trima(i,j+1)==0
                Trimint(i,j)=0;
            end
        end
    end
end
Trima=Trimint;

for i =1:(S(1)-1)
    for j= 1:S(2)
        if Trimb(i,j)==0 && Trimb(i+1,j)>0
            Trimintb(i+1,j)=0;
        else if Trimb(i,j)>0 && Trimb(i+1,j)==0
                Trimintb(i,j)=0;
            end
        end
    end
end
Trimb=Trimintb;
end

Trimed=Trima.*Trimb;



             
