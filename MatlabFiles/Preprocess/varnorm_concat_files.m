clear



for j=1:size(files,2)
clear D Data S R Rvarnorm R2
D=MRIread(files{j});
Data=D.vol;
S=size(Data);

parfor i=1:S(4)
R(:,i)=reshape(Data(:,:,:,i),[1,S(1)*S(2)*S(3)]);
end

%add floats or double here?
 P=prctile(reshape(R,size(R,1)*size(R,2),1),90);
 R1=double(R).*(10000/P);

% demean
S2=size(R1);
parfor i=1:S2(1)
    R2(i,:)=R1(i,:)-mean1(R(i,:));
    %R2(i,:)=medfilt3(R2(i,:),5);
end

%need timeXspace
R=R2';

% % %simple std norm (choose between this and Smith lines below)
stddevs=max(std(R),0.001);  
R=R./repmat(stddevs,size(R,1),1);  % var-norm

%pseudo time series for fsl

for i=1:S(4)
    Rvarnorm(:,:,:,i)=reshape(R(i,:),[S(1),S(2),S(3)]);
end


if exist('Dfinal','var') == 1
    Dfinal=cat(4,Dfinal,Rvarnorm);
else
    Dfinal=Rvarnorm;
end

end

% save
out4.vol = Dfinal;
err = MRIwrite(out4,strcat(files{1}(1:size(files{1},2)-4),'seriesvarnorm.nii'));