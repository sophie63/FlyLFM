
clear all
prompt = 'What is the number of planes?';
Nz = input(prompt)
% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');

file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);

imgdata=reshape(D,[S1(1),S1(2)/Nz,Nz,S1(3)]);

% t = Tiff(strcat(file(1:size(file,2)-4),'2D.tiff'),'w');
% tagstruct.ImageLength = size(imgdata,1)
% tagstruct.ImageWidth = size(imgdata,2)
% tagstruct.Photometric = Tiff.Photometric.RGB
% tagstruct.BitsPerSample = 32
% tagstruct.SamplesPerPixel = S1(4)
% tagstruct.RowsPerStrip = 16
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky
% tagstruct.Software = 'MATLAB'
% t.setTag(tagstruct)
% 
% t.write(imgdata);
% t.close();
out.vol = imgdata;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'3D.nii'));
