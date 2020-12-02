% Downloads neuron images and transforms them to dramatically reduce the size before saving
%VFB_00014201
% need also to get 0000, ....
for k=0:9999
for i=0:9999
     try
 
    url=strcat('https://v2.virtualflybrain.org/data/VFB/i/',sprintf('%04d',k),'/',sprintf('%04d',i),'/VFB_00101567/volume.nrrd')
    outfilename = websave('volume.nrrd',url)
    [X, meta] = nrrdread(outfilename);
    Ds=BoxavFunction3D(X);
    S=size(Ds);
    file2=strcat("Z:\GrunwaldKadow_group\Sophie\WalkProject\VFBData\VFB_",sprintf('%04d',k),sprintf('%04d',i),".tif");  
    imwrite(Ds(:,:,1),file2)
    for j=2:S(3)
    imwrite(Ds(:,:,j),file2,'WriteMode','append')
    end
    delete(outfilename)
      catch 
    % Nothing to do
  end
end
end

