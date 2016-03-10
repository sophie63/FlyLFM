
for i=1:2

  outputFileName = sprintf('img_stack_%d.tif',i)
  file=strcat(folder, outputFileName)
   Duin=uint16(D(:,:,:,i)*32768/2+32768/2);
   Mint=uint16(M);
   Duint3=Duin.*Mint;
for K=1:length(D(1, 1, :,i))
    
   imwrite(Duint3(:, :, K), file, 'WriteMode', 'append');
end  
    
end

    
    