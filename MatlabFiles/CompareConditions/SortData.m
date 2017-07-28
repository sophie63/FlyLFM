%choose the reconstructed 4D nifti file
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);


%Average walk and not walk
Walk=Left+Right+Straight;


%plot and choose a threshold
plot(Walk)

%Rest: if Walk is below threshold
Rest=Walk;
Rest(Walk<50)=50;
Rest(Walk>50)=0;

DataWalk=Data(:,:,:,logical(abs(1-Rest/50)));
DataRest=Data(:,:,:,logical(Rest/50));

out.vol=DataWalk;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'WalkFull.nii'));

out.vol=DataRest;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestFull.nii'));


