[FileName,PathName] = uigetfile('*.nii','Select the registered template file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Masks=D.vol;
Sm=size(Masks);

load('/home/sophie/RegionCorrespondanceForMasks.mat')

Masks2=zeros(Sm(1),Sm(2),Sm(3),87);

parfor j=1:87
    Masks2(:,:,:,j)=(Masks==j);
    j
end

Masks3=zeros(Sm(1),Sm(2),Sm(3),12);

for i=1:12
    for j=1:74
        if LargeRegion(j)==i
            Masks3(:,:,:,i)=Masks2(:,:,:,SmallRegion(j))+Masks3(:,:,:,i);
        end
    end
end
            

out.vol=Masks3;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'separateLarge.nii'));