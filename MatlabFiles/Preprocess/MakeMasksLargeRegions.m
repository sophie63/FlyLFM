MasksLargeRegion=zeros(S(1),S(2),S(3),max(MaskId));
for i=2:86
  for j=1:max(MaskId)
      if MaskId(i+1)==j
        MasksLargeRegion(:,:,:,j)=MasksLargeRegion(:,:,:,j)+Data(:,:,:,i);
      end
  end
end

out.vol=MasksLargeRegion;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'LargeRegions.nii'));

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data2=D.vol;
S=size(Data2);

TSregion=zeros(size(MasksLargeRegion,4),S(4));
for j=1:S(4)
for i=1:size(MasksLargeRegion,4)
    TsRegion(i,j)=sum(sum(sum(MasksLargeRegion(:,:,:,i).*Data2(:,:,:,j))))/sum(sum(sum(MasksLargeRegion(:,:,:,i))));
end
end

figure 
hold on
TsRegionNorm=zeros(size(MasksLargeRegion,4),S(4));
for i=1:size(MasksLargeRegion,4)
    TsRegionNorm(i,:)=TsRegion(i,:);
    plot(1.5*TsRegionNorm(i,:)'+i)
end


