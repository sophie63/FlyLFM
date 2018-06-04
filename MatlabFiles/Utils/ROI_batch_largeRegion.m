clear all
fmask = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_mask_largeRegion_lists.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_dff_lists.txt');
Dlist = strsplit(fdata);
fname = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_name_lists.txt');
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
            TS(j,i)=sum(sum(sum(M2(:,:,:,j).*squeeze(Data(:,:,:,i)))));
        end
    end

    plotTS
    saveas(f, strcat(Nlist{idx},'RegionTS.png'));
    savefig(strcat(Nlist{idx},'RegionTS.fig'));
    save(strcat(Nlist{idx},'TS.mat'),'TS');
    clear TS
    clear Mask
    clear Data
    clear M2
    close all
end