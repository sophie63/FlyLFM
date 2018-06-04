
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


[FileName,PathName] = uigetfile('*.nii','Select the masks Nifti file','/media/sophie/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

Dtrim=Data;

for i=1:S(4)
Z0=1:size(Data,3);
Zinit=((Z0-z1)*dz+dz/2);
%z_psf=abs(Zinit)*0.239+5.46;
%x_psf=abs(Zinit)*0.075+3.3;
% Average the stack layers when the sampling is below psf half width

for j=1:S(3)
nxy=int8((abs(Zinit(j))*0.075+3.3)/(2*dx))
Dtrim(:,:,j,i)=Trim2D(Data(:,:,j,i),nxy);
end

end

out.vol=Dtrim;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'fulltrimmed.nii'));

