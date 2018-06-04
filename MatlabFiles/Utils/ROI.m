%clear all
%open mask
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/media/sophie/008C0665790F0763/ComboPanNeuronalGCaMP6');
file=strcat(PathName,FileName)
M=MRIread(file);
Mask=M.vol;

%Mask=1-Mask;%use that line if mask from segmentation editor
%Md=double(M);

%open data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file','/media/sophie/008C0665790F0763/ComboPanNeuronalGCaMP6');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;

S=size(Data);
M2=Mask./max((max(max(max(Mask)))));
SM2=size(M2);
Data(isnan(Data))=0;

for i=1:S(4)
    for j=1:SM2(4)
        TS(i,j)=sum(sum(sum(M2(:,:,:,j).*squeeze(Data(:,:,:,i)))))/sum(sum(sum(M2(:,:,:,j))));
    end
end


