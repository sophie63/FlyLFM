%parameters to fill up
dz=6;
z1=16;


[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);

for i=1:S(4)
Z0=1:size(Data,3);
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;

% Average the stack layers when the sampling is below psf half width
j=1;
k=1;
while k<=(size(Data,3))
    if (z_psf(k)>dz)
        nz=int8(z_psf(k)/(dz));
        Dpsf2(:,:,j,i)=mean(Data(:,:,k:min((k+nz),size(Data,3)),i),3);
        Znew(j)=Zinit(k)+z_psf(k)/(2*dz);
        j=j+1;
        k=k+nz+1;
    else
        Dpsf2(:,:,j,i)=Data(:,:,k,i);
        Znew(j)=Zinit(k);
        j=j+1;
        k=k+1;
    end
end

i
end

out.vol=Dpsf2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'psf.nii'));
%DataNyq=cell(S(3),1);
%average in x , y if sample is below psf half heigth / 2
%x_psf=abs(Znew)*0.075+3.3;
% for j=1:S(3)
%     nx(j)=int8(x_psf(j)/(2*6/ss));
%     if (nx(j)>1)
%         clear D D1 D2 D3
%         for k=1:S(1)
%             D(1:S(2))=smooth(Dataz(k,:,j),nx(j));
%             D1(k,1:S(2))=downsample(D,nx(j));
%         end
%         for k=1:S(2)/nx(j)
%             D2=smooth(D1(:,k,j),nx(j));
%             D3(1:S(1)/nx(j),1:S(2)/nx(j))=downsample(D2,nx(j));
%         end
%         DataNyq{j,1}=D3;
%     else
%         DataNyq{j,1}=Dataz(:,:,j);
%     end
%     j
% end
