
S=size(TSo);
w=zeros(1,S(2));
%Zo=Left-Right+0.001;
Zo=Right-Left+0.001;
z=sign(Zo);

D=TSo;

for k=1:1000
for i=1:S(1)
    y=w*D(i,:)';
    if sign(y)~=z(i)
        w=w+z(i)*D(i,:);
    end
end
E(k)=sum(z'-w*D');
end

figure
plot(E)



