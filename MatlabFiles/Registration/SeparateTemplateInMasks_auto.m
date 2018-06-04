fid = fopen('/media/sophie/db554c18-e3eb-41e2-afad-7de1c92bf4a5/paper_dataset_lists.txt');
file = fgetl(fid);
while ischar(file)
    file
    D=MRIread(file);
    Masks=D.vol;
    Sm=size(Masks);


    Masks2=zeros(Sm(1),Sm(2),Sm(3),87);

    parfor j=1:87
        Masks2(:,:,:,j)=(Masks==j);
        j
    end


    out.vol=Masks2;
    err = MRIwrite(out,strcat(file(1:size(file,2)-4),'separate.nii'));
    file = fgetl(fid);
end