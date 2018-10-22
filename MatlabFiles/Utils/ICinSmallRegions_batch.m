fmask = fileread('/media/test5/ComponentsList.txt');
Mlist = strsplit(fmask);
fdata = fileread('/media/test5/MasksFreeWalk.txt');
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
    
    D2=Data./(max(max(max(max(Data)))));
    SM2=size(Mask);

    for i=1:S(4)
        for j=1:SM2(4)
            Av(i,j)=sum(sum(sum(Mask(:,:,:,j).*squeeze(D2(:,:,:,i)))))/sum(sum(sum(Mask(:,:,:,j))));
        end
    end

    save(strcat(file(1:size(file,2)-4),'ICsmallRegions.mat'),'Av');
    
    clear Av
    clear Mask
    clear Data
    clear M2
    close all
end