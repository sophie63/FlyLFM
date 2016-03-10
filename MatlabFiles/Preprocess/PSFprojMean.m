function [DataNyq,Znew] = PSFprojMean(Data, dz, z1)

Z0=1:size(Data,3);
Zinit=((Z0-z1)*dz+dz/2);
z_psf=abs(Zinit)*0.239+5.46;

% Average the stack layers when the sampling is below psf half width
j=1;
i=1;
while i<=(size(Data,3))
    if (z_psf(i)>dz)
        nz=int8(z_psf(i)/(dz));
        Dataz(:,:,j)=mean(Data(:,:,i:min((i+nz),size(Data,3))),3);
        Znew(j)=Zinit(i)+z_psf(i)/(2*dz);
        j=j+1;
        i=i+nz+1;
    else
        Dataz(:,:,j)=Data(:,:,i);
        Znew(j)=Zinit(i);
        j=j+1;
        i=i+1;
    end
end

%average in x , y if sample is below psf half heigth / 2
x_psf=abs(Znew)*0.075+3.3;

S=size(Dataz);
DataNyq=Dataz;
%DataNyq=cell(S(3),1);

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
