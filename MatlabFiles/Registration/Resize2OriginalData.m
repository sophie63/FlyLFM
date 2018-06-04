[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=squeeze(D.vol);
S=size(Data);

Datar=imresize3d(Data,[],[83 165 71],'linear');

out.vol=Datar;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'r.nii'));