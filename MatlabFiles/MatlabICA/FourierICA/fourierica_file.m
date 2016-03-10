clear

%open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=B.vol;
S1=size(D)


parfor i=1:S1(4)
R(:,i)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
end

clear D
clear B

[S_FT,A_orig,W_orig]=fourierica(R,150,200,0.01,90);

 for i=1:150
Maps(:,:,:,i)=reshape(W_orig(i,:),[S1(1),S1(2),S1(3)]);
end
size(Maps)



out.vol=Maps;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'FourierMaps.nii'));

