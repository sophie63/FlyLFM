%make average of the time series

%javaaddpath '/home/sophie/Matlab/java/mij.jar'
javaaddpath '/usr/local/MATLAB/R2014b/java/mij/mij.jar'
javaaddpath '/home/sophie/Fiji.app/jars/ij-1.48q.jar'
MIJI;
%drag and drop data
Data=MIJ.getCurrentImage;

% use the subset of LPUs that are clearly visible is data (4landmarks.tif),
% record and cropp in FIJI to have grosso modo the same volume as the data

% I1 is the template that has been cropped to look like the data I2
I1=double(I1); I2=double(I2);

% First resize volume I1 to match size of volume I2
I1r=imresize3d(I1,[],size(I2),'linear');

% resize all the LPUs infos so that they match the Data
S=size(CLcropped)
for i=1:S(4)
    CLr(:,:,:,i)=imresize3d(CLcropped(:,:,:,i),[],size(I2),'linear');
end

% save I1 and open in FIJI to get the landmarks
MIJ.createImage(I2)
MIJ.createImage(I1)

 % Me is the landmarks for my data, Chiang for the template data
 % need to inverse x and y
 
 [O_transMC2,Spacing,Xreg]=point_registration([98 199 93],Me,Chiang);
 IregMC2=bspline_transform(O_transMC2,I1r,Spacing,3);
 MIJ.createImage(IregMC2)
 MIJ.createImage(I2)
 
for i=1:S(4)
    CLreg(:,:,:,i)=bspline_transform(O_transMC2,CLr(:,:,:,i),Spacing,3);
end


 % Start again to refine the landmarks 

 [O_transMC2b,Spacingb,Xregb]=point_registration([98 199 93],Meb,Chiangb);
 IregMC2b=bspline_transform(O_transMC2b,I1r,Spacingb,3);
 MIJ.createImage(IregMC2b)
 MIJ.createImage(I2)
 
for i=1:S(4)
    CLregb(:,:,:,i)=bspline_transform(O_transMC2b,CLreg(:,:,:,i),Spacingb,3);
end


%Binarize
CLregbin=CLregb;
CLregbin(CLreg>0.0001)=1;

%import time series
 DataTS=MIJ.getCurrentImage;
DTS=reshape(DataTS,98,199,93,1159);
S2=size(DTS)
x=1:S2(4);

%play with the initial values of the fit for bleaching to have nice fits
for i=1:S(4)
%         for (j=1:S2(4))
%     D(:,:,:,j)=DTS(:,:,:,j).*CLregbin(:,:,:,i);
%     Mean_Data_in_LPUs(j,i)=mean(mean(mean(D(:,:,:,j))));
%     end   
     A=fit(x',Mean_Data_in_LPUs(:,i),'a*(1+b*exp(c*x)-0.29*exp(0.00012*x))','StartPoint',[0.0001,0.045,-0.03],'Lower',[0.00001,-5,-2],'Upper',[100,2,1])
%     
    figure(i)
    plot(A,x',Mean_Data_in_LPUs(:,i))
    B=squeeze(A(x));
    U_LPUs(:,i)=Mean_Data_in_LPUs(:,i)-B(:);
end

LPUactivity=-U_LPUs;
for i=1:S(4)
    figure(i)
    plot(LPUactivity(:,i))
end



