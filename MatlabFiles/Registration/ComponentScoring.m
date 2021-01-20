St=[242 113 34]

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=double(D.vol);
Datab = permute(Data,[2 1 3 4]);
S=size(Datab);

%Open the goodICs .mat
[FileName,PathName] = uigetfile('*.mat','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
load(file)

%Open the correspondance points .mat
[FileName,PathName] = uigetfile('*.mat','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
load(file)

%B-spline
%[O_transMC2,Spacing,Xreg]=point_registration([113 242 34],TempPoints(:,idx),DataPoints(:,idx));
%IregMC2=bspline_transform(O_transMC2,Data(:,:,:,16),Spacing,3);


%Thin plate
for x=1:St(1)
for y=1:St(2)
for z=1:St(3)
[wobject] = TPS3D(TempPoints,DataPoints,[x y z]);
zobject(:,z)=wobject;
end
for n=1:size(GoodIC,2)
NewVals(x,y,:,n)=interp3(Data(:,:,:,GoodIC(n)),zobject(1,:),zobject(2,:),zobject(3,:));
end

end
x
end

for n=1:size(GoodIC,2)
NewValsc(:,:,:,n)=BoxavFunction3Dby2(NewVals(:,:,:,n));
end

%out.vol=IregMC2;
NewValsb=permute(NewValsc,[2 1 3 4]);
out.vol=NewValsb;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'Tempreg.nii'));



