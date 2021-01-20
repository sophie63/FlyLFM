clear all

Folders=dir("/media/NAS/Sophie/RawData/ss1");

for l=3:size(Folders)   

files=dir(strcat("/media/NAS/Sophie/RawData/ss1/",Folders(l).name));
Folders(l).name

for j=3:size(files,1)
    FileName=strcat(Folders(l).name,'-',num2str((j-2),'%05.f'),'.tif')
tiff_info = imfinfo(strcat("/media/NAS/Sophie/RawData/ss1/",Folders(l).name,"/",FileName)); % return tiff structure, one element per image
tiff_stack = imread(strcat("/media/NAS/Sophie/RawData/ss1/",Folders(l).name,"/",FileName), 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(strcat("/media/NAS/Sophie/RawData/ss1/",Folders(l).name,"/",FileName), ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
Stack(:,:,:,j)=tiff_stack;

clear tiff_stack
end

idx=[2 1 3];

out.vol=Stack;
err = MRIwrite(out,strcat('/media/NAS/Sophie/RawData/toPreprocess/',Folders(l).name,'.nii'));
clear Stack out
end

