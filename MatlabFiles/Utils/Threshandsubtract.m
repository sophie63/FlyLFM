
% Open Data
[FileName,PathName] = uigetfile('*.nii','Select the Nifti file');
file=strcat(PathName,FileName)
B=MRIread(file);
D=double(B.vol);
clear B
S1=size(D);


parfor i=1:S1(4)
R(i,:)=reshape(D(:,:,:,i),[1,S1(1)*S1(2)*S1(3)]);
S=std(R(i,:));
R2=R(i,:)-3*S;
R2(R2<0)=0;
D2(:,:,:,i)=reshape(R2,[S1(1),S1(2),S1(3)]);
end
out.vol = D2;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'thresh3stdmin3std.nii'));