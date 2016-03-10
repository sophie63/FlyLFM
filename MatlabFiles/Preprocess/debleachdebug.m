for (i=1:S(1))
    i
    for (j=20:S(2))
        for(k=10:20:2)  
        D=squeeze(Data(i,j,k,:));
     %   A=fit(t',D,'exp2','StartPoint',[1,-0.001,0.1,-0.001],'Lower',[0.0001,-0.00001,0.00001,-0.00001]);\
      %Use values from fir of the average of all points
        A=fit(t',D,'a*(1+b*exp(-c*x)+d*exp(-e*x))','StartPoint',[mean(mean(mean(D))),0.0556,0.01377,0.1629,0.001028],'Lower',[mean(mean(mean(D)))/4,0.01,0.002,0.01,0],'Upper',[mean(mean(mean(D)))*2,0.6,0.1,0.6,0.002]);
        B=squeeze(A(t));

        Unbleached_data(i,j,k,:)=D-B;
        if(j==25)
            figure(i)
            plot(A,t',D)
        end
        end
        
    end
end
