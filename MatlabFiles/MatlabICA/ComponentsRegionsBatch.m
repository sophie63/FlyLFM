fmask = fileread('/media/test15/JFRCpsfList.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/test15/ComponentsDataList.txt');
Dlist = strsplit(fdata);

for idx=1:length(Mlist)
  
    M=MRIread(Mlist{idx});
    Mask=M.vol;

    file=Dlist{idx}
    file2=strcat(file(1:size(file,2)-4),'thresh3std.nii');
    D=MRIread(file2);
    Data=D.vol;
    
    S=size(Data);
    M2=Mask./(max(max(max(max(Mask)))));
    SM2=size(M2);

    for i=1:S(4)
        for j=1:SM2(4)
            Sm(i,j)=sum(sum(sum(M2(:,:,:,j).*squeeze(Data(:,:,:,i)))))/sum(sum(sum(Data(:,:,:,i))));
        end
    end

    save(strcat(file(1:size(file,2)-4),'sumICSegions.mat'),'Sm');
    
    clear Sm
    clear Mask
    clear Data
    clear M2
    close all
    clear M
    clear D
end