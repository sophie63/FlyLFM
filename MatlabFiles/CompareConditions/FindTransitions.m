j=1

for i=20:size(Walk,2)-20
    if (all(Walk(i-19:i)<10) & all(Walk(i+1:i+20)>10))
        TransitionUp(j)=i;
        j=j+1;
    end
end

j=1
for i=20:size(Walk,2)-20
    if (all(Walk(i-19:i)>10) & all(Walk(i+1:i+20)<10))
        TransitionDown(j)=i;
        j=j+1;
    end
end


DataTrUp=zeros(S(1),S(2),S(3),100);
DataTrUpb=zeros(S(1),S(2),S(3),50);

for i=1:size(TransitionUp)
    DataTrUp=DataTrUp+Data(:,:,:,(TransitionUp(i)-50):(TransitionUp(i)+49));
    DataTrUpb=DataTrUpb+Data(:,:,:,(TransitionUp(i)-49):TransitionUp(i));
end


out.vol=DataTrUp;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'TransitionUp.nii'));

out.vol=DataTrUpb;
err = MRIwrite(out,strcat(file(1:size(file,2)-4),'BeforeTransitionUp.nii'));
