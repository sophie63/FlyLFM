%parameters to fill up
dz=6;
z1=16;


[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);


Z0=1:size(Data,3);
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;

% Average the stack layers when the sampling is below psf half width
j=1;
k=1;
while k<=(size(Data,3))
    if (z_psf(k)>dz)
        nz=int8(z_psf(k)/(dz));
        Dpsf2(:,:,j)=mean(Data(:,:,k:min((k+nz),size(Data,3))),3);
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j)=Data(:,:,k);
        Znew(j)=Zinit(k);
        j=j+1;
        k=k+1;
    end
end

out.vol=Dpsf2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'psf.nii'));

