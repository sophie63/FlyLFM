
clear all

fname = fileread('/media/test15/ComponentsDataList.txt');
Nlist = strsplit(fname);

for idx=1:length(Nlist)

    Data=MRIread(Nlist{idx});
    D=Data.vol;
    
    S1=size(D);

    parfor i=1:S1(4)
    R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
    S=std(R(i,:));
    R2=R(i,:);
    R2(abs(R2)<(3*S))=0;
    D2(:,:,:,i)=reshape(R2,[S1(1),S1(2),S1(3)]);
    end
    out.vol = D2;
    file=Nlist{idx}
    err = MRIwrite(out,strcat(file(1:size(file,2)-4),'thresh3std.nii'));
   
    clear D
    clear R
    clear S
    clear R2
    clear D2
    clear out
end
