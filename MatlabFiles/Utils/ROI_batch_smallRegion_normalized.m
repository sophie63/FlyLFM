clear all
fmask = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_mask_smallRegion_listsV2.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_dff_listsV2.txt');
Dlist = strsplit(fdata);
fname = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_name_listsV2.txt');
Nlist = strsplit(fname);

for idx=1:length(Mlist)
   
    M=MRIread(Mlist{idx});
    Mask=M.vol;
    Mlist{idx}
    
    D=MRIread(Dlist{idx});
    Data=D.vol;
    
    S=size(Data);
    M2=Mask./(max(max(max(max(Mask)))));
    SM2=size(M2);

    for i=1:S(4)
        for j=1:SM2(4)
            TS(j,i)=sum(sum(sum(M2(:,:,:,j).*squeeze(Data(:,:,:,i)))))/sum(sum(sum(M2(:,:,:,j))));
        end
    end

    save(strcat(Nlist{idx},'smallRegions_TS_normalized.mat'),'TS');
    clear TS
    clear Mask
    clear Data
    clear M2
    close all
end