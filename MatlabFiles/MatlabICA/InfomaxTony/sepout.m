% SEPOUT - put whatever textual output report you want here.
%  Called after each pass through the data.
%  If your data is real, not artificially mixed, you will need
%  to comment out line 4, since you have no idea what the matrix 'a' is.
% 
[change,olddelta,angle]=wchange(oldw,w,olddelta); 
oldw=w;
fprintf('****sweep=%d, change=%.4f angle=%.1f deg., [N%d,M%d,P%d,B%d,L%.5f] \n',...
   sweep,change,180*angle/pi,N,M,P,B,L);
w*wz*a     %should be a permutation matrix for artif. mixed data.
