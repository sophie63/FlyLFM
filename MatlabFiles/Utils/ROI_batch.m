fmask = fileread('/media/test5');
Mlist = strsplit(fmask);
fdata = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_dff_lists.txt');
Dlist = strsplit(fdata);
fname = fileread('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_name_lists.txt');
Nlist = strsplit(fname);

for idx=1:length(Mlist)
   
    M=MRIread(Mlist{idx});
    Mask=M.vol;

    D=MRIread(Dlist{idx});
    Data=D.vol;
    
    S=size(Data);
    M2=Mask./(max(max(max(max(Mask)))));
    SM2=size(M2);


    M2r=reshape(M2,SM2(1)*SM2(2)*SM2(3),SM2(4));
    Datar=reshape(Data,SM2(1)*SM2(2)*SM2(3),S(4));
    M2rnorm=normc(M2r);

    TS=M2rnorm'*Datar;
    plotTS
    saveas(f, strcat(Nlist{idx},'RegionTS_LR.png'));
    savefig(strcat(Nlist{idx},'RegionTS_LR.fig'));
    save(strcat(Nlist{idx},'TS_LR.mat'),'TS');
    clear TS
    clear Mask M2r M2rnorm
    clear Data Datar
    clear M2
    close all
end