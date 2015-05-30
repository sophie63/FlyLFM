function [z,w,h]=PSF(P,D)

% This function measure the psf half heigt width along x and z
%
% Inputs 
% P: coordinates of the center of the small fluorescent bead
% (determined with imageJ)
% D: Whole image of fluorescent beads embedded in a gel
%
% Output: 
% z : distance of the bead center to the focal plane
% w : lateral half heigt width
% h : axial half heigt width

% example use: 
% P=[100,120,20];
% [z(1),w(1),h(1)]=PSF(P,D)
% repeat for beads at different depth


Dp1=D(P(1)-6:P(1)+6,P(2)-6:P(2)+6,P(3)-10:P(3)+10);
Xp1=squeeze(max(max(Dp1,[],3),[],2));
Xp1=Xp1-(Xp1(1)+Xp1(1:13))/2;
plot(squeeze(Xp1));
title('X');
[M,xm]=max(squeeze(Xp1));
xa=interp1(Xp1(1:xm),[1:xm],M/2);
xb=interp1(Xp1(xm:13),[xm:13],M/2);
w=xb-xa;

Zp1=squeeze(max(max(Dp1,[],1),[],2));
Zp1=Zp1-(Zp1(1)+Zp1(20))/2;
figure
plot(squeeze(Zp1))
[M,z]=max(squeeze(Zp1))
za=interp1(Zp1(1:z),[1:z],M/2);
zb=interp1(Zp1(z:20),[z:20],M/2);
h=zb-za;
z=P(3)-10+z;



