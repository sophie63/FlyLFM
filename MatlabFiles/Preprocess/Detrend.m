function Unbleached_data = Detrend(Data ,Fr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

S=size(Data)
t=1:S(4);
%Unbleached_data=zeros(S);

for (i=1:S(1))
    i
    for (j=1:S(2))
        parfor(k=1:S(3))
        D=squeeze(Data(i,j,k,:));
    
%Use values from fir of the average of all points
        %Dm=mean(D)
        %A=fit(t',D,'a*(1+b*exp(-c*x)+d*exp(-e*x))','StartPoint',[Dm,0.0556,0.007,0.1629,0.0005],'Lower',[Dm/4,0.01,0.001,0.01,0],'Upper',[Dm*2,0.6,0.05,0.6,0.002]);
<<<<<<< HEAD
        B=smooth(D,14/Fr);
||||||| merged common ancestors
        B=smooth(D,7/Fr);
=======
        % Try also median first: very slow so remove the following lines if in a hurry
        %B=medfilt1(D,10/Fr);
        C=smooth(D,10/Fr);
>>>>>>> 391ef8394753bc7dd73a58e5f42c106fc9574953

        Unbleached_data(i,j,k,:)=D-C;
%         if(j==25)
%             figure(i)
%             plot(A,t',D)
%         end
        end

    end
end
