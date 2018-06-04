
% This code trims the anatomical masks to avoid crosscontamination of the
% brain regions due to low resolution

clear all
%parameters to fill up

prompt = 'What is the position of the focal plane in the stack?';
z1 = input(prompt)

prompt = 'What is the distance betweeen stacks?';
dz= input(prompt)

prompt = 'What is the xy resolution?';
dx1= input(prompt)

prompt = 'What is the magnification?';
Mag= input(prompt)

dx=dx1

[FileName,PathName] = uigetfile('*.nii','Select the average of the data Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
data=D.vol;
S=size(data);

[FileName,PathName] = uigetfile('*.nii','Select the registered template file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Masks=D.vol;
Sm=size(Masks);


Masks2=zeros(Sm(1),Sm(2),Sm(3),87);
MasksResized=zeros(S(1),S(2),S(3),87);
for j=1:87
    Masks2(:,:,:,j)=(Masks==j);
    parfor i=1:Sm(3)
        MasksResized(:,:,i,j)=imresize(Masks2(:,:,i,j),[S(1),S(2)]);
    end
end

Data=MasksResized;
S=size(MasksResized);

for i=1:S(4)
Z0=1:size(Data,3);
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;
x_psf=abs(Zinit)*0.075+3.3;
% Average the stack layers when the sampling is below psf half width
j=1;
k=1;
while k<=(size(Data,3))
    if (z_psf(k)>dz)
        nz=int8(z_psf(k)/(dz));
        Dpsf2(:,:,j,i)=Data(:,:,min((k+int8(nz/2)),size(Data,3)),i);
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        nxy=int8((abs(Znew(j))*0.075+3.3)/(2*dx))
        Dtrim(:,:,j,i)=Trim2D(Dpsf2(:,:,j,i),nxy);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j,i)=Data(:,:,k,i);
        Znew(j)=Zinit(k);
        nxy=int8((abs(Znew(j))*0.075+3.3)/(2*dx));
        Dtrim(:,:,j,i)=Trim2D(Dpsf2(:,:,j,i),nxy);
        j=j+1;
        k=k+1;
    end
end

end

out.vol=Dtrim;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'fullpsftrimmed.nii'));