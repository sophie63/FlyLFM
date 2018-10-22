fmask = fileread('/media/test15/ComponentsDataList2.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/test15/DataForROI2.txt');
Dlist = strsplit(fdata);

for idx=1:length(Mlist)
   
    file=Mlist{idx}
    file2=strcat(file(1:size(file,2)-4),'thresh3std.nii');
    M=MRIread(file2);
    Mask=M.vol;
    Dlist{idx}
    D=MRIread(Dlist{idx});
    Data=D.vol;
    
    S=size(Data);
    M2=Mask./(max(max(max(max(Mask)))));
    SM2=size(M2);

    for i=1:S(4)
        for j=1:SM2(4)
            TS(i,j)=sum(sum(sum(M2(:,:,:,j).*squeeze(Data(:,:,:,i)))))/sum(sum(sum(M2(:,:,:,j))));
        end
    end

    save(strcat(file(1:size(file,2)-4),'TSROIThresh3.mat'),'TS');
    
    clear TS
    clear Mask
    clear Data
    clear M2
    close all
end