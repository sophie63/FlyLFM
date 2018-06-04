% This code trims the anatomical masks to avoid crosscontamination of the
% brain regions due to low resolution

clear all
%parameters to fill up

prompt = 'What is the position of the focal plane in the stack?';
z1 = input(prompt)

prompt = 'What is the distance betweeen stacks?';
dz= input(prompt)

prompt = 'What is the xy resolution?';
dx= input(prompt)


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

S=size(Masks3);

for i=1:S(4)
Z0=1:size(Masks3,3);
Zinit=((Z0-z1)*dz+dz/2);
%z_psf=abs(Zinit)*0.239+5.46;
%x_psf=abs(Zinit)*0.075+3.3;
% Average the stack layers when the sampling is below psf half width

for j=1:S(3)
nxy=int8((abs(Zinit(j))*0.075+3.3)/(2*dx));
Dtrim(:,:,j,i)=Trim2D(Masks3(:,:,j,i),nxy);
end

end

out.vol=Dtrim;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'separateLargetrimmed.nii'));