fmask = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_mask_largeRegion_lists.txt');
Mlist = strsplit(fmask);

fname = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_name_lists.txt');
Nlist = strsplit(fname);

for idx=1:length(Mlist)
   
    M=MRIread(Mlist{idx});
    Mask=M.vol;
    Mlist{idx}

    M2=Mask./(max(max(max(max(Mask)))));
    SM2=size(M2);


    for j=1:SM2(4)
            TSnorm(j)=sum(sum(sum(M2(:,:,:,j))));
    end


    save(strcat(Nlist{idx},'TSnorm.mat'),'TSnorm');
    clear TSnorm
    clear Mask
    clear M2

end