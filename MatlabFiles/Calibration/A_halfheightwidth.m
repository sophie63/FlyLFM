function [z,w,h]=PSF(P)

Dp1=D(P(1)-5:P(1)+5,P(2)-5:P(2)+5,P(3)-5:P(3)+5);
Xp1=max(max(Dp1,[],3),[],2);
[M,xm]=max(Xp1);
xa=interp1(Xp1(1:xm),[1:xm],M/2);
xb=interp1(Xp1(xm:11),[xm:11],M/2);
w=xb-xa


