[FileName,PathName] = uigetfile('*.nii','Select the registered template file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Masks=D.vol;
Sm=size(Masks);


Masks2=zeros(Sm(1),Sm(2),Sm(3));
M=squeeze(max(max(max(max(Masks)))));

parfor j=1:87
    Masks2=Masks2+j*Masks(:,:,:,j)/M;
    j
end


out.vol=Masks2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'together.nii'));