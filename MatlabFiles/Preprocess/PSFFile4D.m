% This code trims the anatomical masks to avoid crosscontamination of the
% brain regions due to low resolution

clear all

%parameters to fill up
dz=12;
prompt = 'What is the position of the focal plane in the stack?';
z1 = input(prompt)

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/media/sophie/1554f3a2-94a1-4cf8-ae79-8a6fef1b5f7e/FreeBehaviorPanNeuronalGCaMP6/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);
clear D


Z0=1:size(Data,3);
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;
% Average the stack layers when the sampling is below psf half width

for l=1:S(4)
j=1;
k=1;
while k<=(size(Data,3))
    if (z_psf(k)>dz)
        nz=int8(z_psf(k)/(dz));
        Dpsf2(:,:,j,l)=Data(:,:,min((k+int8(nz/2)),size(Data,3)),l);
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j,l)=Data(:,:,k,l);
        Znew(j)=Zinit(k);
        j=j+1;
        k=k+1;
    end
end
end


out.vol=Dpsf2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'psf.nii'));
