 B=MRIread('~/Desktop/779ss1regUc_250.ica/filtered_func_data.ica/779ss1regUcmelodic_IC.nii');
 D=B.vol;
folder = uigetdir
S=size(D)

A=max(max(max(max(D))))


for i=1:S(4)

  outputFileName = sprintf('img_stack_%03d.tif',i)
   Din=int16(D(:,:,:,i)*65536/A);
   
  t = Tiff(outputFileName,'a');
   tagstruct.ImageLength = size(Din,1)
   tagstruct.ImageWidth = size(Din,2)
tagstruct.Photometric = Tiff.Photometric.MinIsBlack
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky
tagstruct.BitsPerSample   = 32;
tagstruct.Software = 'MATLAB'
t.setTag(tagstruct)
 
for K=1:length(D(1, 1, :,i))
    t.write(Din(:,:,K));
end  
    
end

    
    