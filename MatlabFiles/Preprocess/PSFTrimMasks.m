<<<<<<< HEAD
% This code trims the anatomical masks to avoid crosscontamination of the
% brain regions due to low resolution

%parameters to fill up
dz=6;
dx=3;
z1=18;

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

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
        Dpsf2(:,:,j,i)=Data(:,:,k+int8(nz/2),i);
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        nxy=int8((abs(Znew(j))*0.075+3.3)/(2*dx))
        Dtrim(:,:,j,i)=Trim2D(Dpsf2(:,:,j,i),nxy);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j,i)=Data(:,:,k,i);
        Znew(j)=Zinit(k);
        nxy=int8((abs(Znew(j))*0.075+3.3)/(2*dx))
        Dtrim(:,:,j,i)=Trim2D(Dpsf2(:,:,j,i),nxy);
        j=j+1;
        k=k+1;
    end
end

end

out.vol=Dtrim;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'fullpsftrimmed.nii'));
||||||| merged common ancestors
=======
% This code trims the anatomical masks to avoid crosscontamination of the
% brain regions due to low resolution

clear all
%parameters to fill up
dz=3;
dx=3;
z1=32;

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

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
        nxy=int8((abs(Znew(j))*0.075+3.3)/(2*dx))
        Dtrim(:,:,j,i)=Trim2D(Dpsf2(:,:,j,i),nxy);
        j=j+1;
        k=k+1;
    end
end

end

out.vol=Dtrim;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'fullpsftrimmed.nii'));
>>>>>>> 57759d29a73933bc92c6aaedf808d46b8fabe631
