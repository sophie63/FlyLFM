function Unbleached_data = dFF(Data,Npoint)
%This function removes a smooth version of the time series (to correct for
%bleaching for example)

S=size(Data)
t=1:S(4);

for (i=1:S(1))
    i
      for(k=1:S(3))
           %parfor (j=1:S(2))
           for (j=1:S(2))
        D=squeeze(Data(i,j,k,:));
        C=smooth(D,Npoint);
        Cm=mean(C);
        Unbleached_data(i,j,k,:)=(D-C)/max([Cm,1]);
        end

    end
end
