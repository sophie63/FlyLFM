folder = uigetdir

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=B.vol;

S=size(D)

for i=1:S(4)

  outputFileName = sprintf('img_stack_%03d.tif',i)
  file=strcat(folder,'/', outputFileName)
   %Duint32=uint32(D(:,:,:,i));
for K=1:length(D(1, 1, :,i))
   imwrite(D(:, :, K,i), file, 'WriteMode', 'append');
end  
    
end

    
    