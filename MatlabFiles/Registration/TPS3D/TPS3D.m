%======================================================
% 3D Thin Plate Spline by Yang Yang (05.05.2012)
% 
%======================================================
function [wobject] = TPS3D(points, ctrlpoints,object)

%======================================================
% Calculate Parameters 
%======================================================
npnts = size(points,1);
K = zeros(npnts, npnts);
for rr = 1:npnts
    for cc = 1:npnts
        K(rr,cc) = sum( (points(rr,:) - points(cc,:)).^2 ); %R^2 
        K(cc,rr) = K(rr,cc);
    end;
end;
%calculate kernel function R
K = max(K,1e-320); 
%K = K.* log(sqrt(K));
K = sqrt(K); %
% Calculate P matrix
P = [ones(npnts,1), points]; %nX4 for 3D
% Calculate L matrix
L = [ [K, P];[P', zeros(4,4)] ]; %zeros(4,4) for 3D
param = pinv(L) * [ctrlpoints; zeros(4,3)]; %zeros(4,3) for 3D

%======================================================
% Calculate new coordinates (x',y',z') for each points 
%====================================================== 
pntsNum=size(object,1); 
K = zeros(pntsNum, npnts);
gx=object(:,1);
gy=object(:,2);
gz=object(:,3);
for nn = 1:npnts
    K(:,nn) = (gx - points(nn,1)).^2 + (gy - points(nn,2) ).^2 + (gz - points(nn,3) ).^2; % R^2
end;
K = max(K,1e-320); 
K = sqrt(K); %|R| for 3D

P = [ones(pntsNum,1), gx, gy, gz];
L = [K, P];
wobject = L * param;  
wobject(:,1)=round(wobject(:,1)*10^3)*10^-3;
wobject(:,2)=round(wobject(:,2)*10^3)*10^-3;
wobject(:,3)=round(wobject(:,3)*10^3)*10^-3;