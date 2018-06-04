% This code trims the anatomical masks to avoid crosscontamination of the
% brain regions due to low resolution

%parameters to fill up

clear all

prompt = 'What is the position of the focal plane in the stack?';
z1 = input(prompt)

prompt = 'What is the distance betweeen stacks?';
dz = input(prompt)

prompt = 'What is the sizeof the initial stack?';
Zfull = input(prompt)

[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
DM=D.vol;
S=size(DM);

Z0=1:Zfull;
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;

% Average the stack layers when the sampling is below psf half width
j=1;
k=1;
while k<=Zfull
    if (z_psf(k)>dz)
        nz=int8(z_psf(k)/(dz));
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        j=j+1;
        k=k+nz+1;
    else
        Znew(j)=Zinit(k);
        j=j+1;
        k=k+1;
    end
end

Z2inta=Zinit(Zinit>min(Znew));
Z2int=Z2inta(Z2inta<max(Znew));

Zbool=(Zinit>min(Znew)).*(Zinit<max(Znew));

Dipsf=zeros(S(1),S(2),size(Zinit,2),S(4));

for i=1:S(1)
    for j=1:S(2)
        for k=(1:S(4))
            Dipsf(i,j,logical(Zbool),k)=interp1(Znew,squeeze(DM(i,j,:,k)),Z2int,[]);
            %Dipsf(i,j,(Zinit>min(Znew)),k)=Dipsf(i,j,Z2int(1),k);
            %Dipsf(i,j,(Zinit<max(Znew)),k)=Dipsf(i,j,Z2int(end),k);
        end
    end
    i
end

out.vol=Dipsf;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'invpsf.nii'));
