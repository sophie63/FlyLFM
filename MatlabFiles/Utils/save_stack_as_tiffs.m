 B=MRIread('713ss1partodor.nii');
 D=B.vol;
folder = uigetdir
S=size(D)

A=max(max(max(max(D))))

for i=1:S(4)

  outputFileName = sprintf('img_stack_%03d.tif',i)
  Din=int16(D(:,:,:,i)*65536/A);
for K=1:length(D(1, 1, :,i))
   imwrite(Din(:, :, K), outputFileName, 'WriteMode', 'append');
end  
    
end

    
    