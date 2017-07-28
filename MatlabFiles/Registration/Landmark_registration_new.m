[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Temp=squeeze(D.vol);
St=size(Temp);

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;

S=size(Data);
idx=[2 1 3];
 [O_transMC2,Spacing,Xreg]=point_registration([256 512 109],JFRCpoints(:,idx),Points945(:,idx));
 IregMC2=bspline_transform(O_transMC2,Data,Spacing,3);
 
out.vol=IregMC2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'Tempreg.nii'));
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

for i=1:S(4)
    CLreg(:,:,:,i)=bspline_transform(O_transMC2,CLr(:,:,:,i),Spacing,3);
end




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



