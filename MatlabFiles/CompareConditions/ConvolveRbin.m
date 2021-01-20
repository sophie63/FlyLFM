clear
file='Z:\n8m\GrunwaldKadow_group\Sophie\WalkProject\Vglut_Gad\Vglut\B167-B169\B167Rbin.mat';
%file='C:\Users\Admin\Downloads\B38Rbin.mat';
load(file)
load('GCaMP6M10ms.mat')
FR=5;
Ratio=100/FR;
Kernel=GCaMP6MKernel(1:Ratio:200);
Rk=conv(R',Kernel);

%dff 4000
Rs=smooth(Rk,4000);
Rd=Rk-Rs';
Rkd=Rd(1:size(R));
plot(Rkd)
save(strcat(file(1:size(file,2)-4),'kd.mat'),'Rkd')