
AvTS=zeros(size(TS,1),601);
for i=1:4
    AvTS=AvTS+TS(:,1+(i-1)*601:i*601);
end
