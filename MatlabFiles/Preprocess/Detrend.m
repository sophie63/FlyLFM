function Unbleached_data = Detrend(Data ,Fr)
%This function removes a smooth version of the time series (to correct for
%bleaching for example)

S=size(Data)
t=1:S(4);

for (i=1:S(1))
    i
    for (j=1:S(2))
        parfor(k=1:S(3))
        D=squeeze(Data(i,j,k,:));
        C=smooth(D,10/Fr);
        Unbleached_data(i,j,k,:)=D-C;
        end

    end
end
