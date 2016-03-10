function Unbleached_data = Debleach( Data )
% This function substracts a double exponential fitted to the data
% The initial values were obtained by trial and error, but seem
% to work for all the GCAMP6 and Arclight data (400hz to 50Hz sampling
% rate)
%
% Input and output: 4D data (x,y,z,t)


S=size(Data)
t=1:S(4);

for (i=1:S(1))
    i
    for (j=1:S(2))
        parfor(k=1:S(3))  
        D=squeeze(Data(i,j,k,:));
     %   A=fit(t',D,'exp2','StartPoint',[1,-0.001,0.1,-0.001],'Lower',[0.0001,-0.00001,0.00001,-0.00001]);
        Dm=mean(D)
        A=fit(t',D,'a*(1+b*exp(-c*x)+d*exp(-e*x))','StartPoint',[Dm,0.0556,0.01377,0.1629,0.001028],'Lower',[Dm/4,0.01,0.002,0.01,0],'Upper',[Dm*2,0.6,0.1,0.6,0.002]);
        B=squeeze(A(t));

        Unbleached_data(i,j,k,:)=D-B;
%         if(j==25)
%             figure(i)
%             plot(A,t',D)
%         end
        end
        
    end
end






