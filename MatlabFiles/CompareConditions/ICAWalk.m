%choose the reconstructed 4D nifti file 
%(Mkf)
[FileName,PathName] = uigetfile('*.nii','Select the Nifti 4D dataset','/home/sophie/Desktop/');
file=strcat(PathName,FileName)
D=MRIread(file);
Data=D.vol;
S=size(Data);


%Average walk and not walk
Walk=Left+Right+Straight;

%plot and choose a threshold
plot(Walk)

%Rest: if Walk is below threshold
Rest=Walk;
Rest(Walk<50)=50;
Rest(Walk>50)=0;

RightId=zeros(1,S(4));
LeftId=zeros(1,S(4));
StraightId=zeros(1,S(4));

for i=1:S(4)
    if ((Left(i)+Right(i)+Straight(i))> 60 & Left(i)>1.5*(Right(i)+Straight(i)))
        LeftId(i)=1;
    else if ((Left(i)+Right(i)+Straight(i))>60 & Right(i)>1.5*(Straight(i)+Left(i)))
        RightId(i)=1;
     else if ((Left(i)+Right(i)+Straight(i))>60 & 1.5*(Left(i)+Right(i))<Straight(i))
        StraightId(i)=1;
         end
        end
    end
end

%out.vol=Data(:,:,:,logical(Rest/50));
%err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RestFull.nii'));

%out.vol=Data(:,:,:,(logical(-(Rest/50-1))));
%err = MRIwrite(out,strcat(file(1:size(file,2)-4),'WalkFull.nii'));

out.vol=Data(:,:,:,logical(LeftId));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'LeftFull.nii'));

out.vol=Data(:,:,:,logical(RightId));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'RightFull.nii'));

out.vol=Data(:,:,:,logical(StraightId));
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'StraightFull.nii'));
  
