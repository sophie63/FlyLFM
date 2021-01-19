function [dN,t]=binspikes(data,Fs,t)
% bin spikes at a specified frequency sampling i.e. sampling rate 1/sampling
% eg: 1ms accuracy use sampling = 1000
% This function mimics the chronux 2.12 function of the same name, with a
% the following bug fixes
%   - handle trials with no spikes: In chronux 2.12, spikeless trials can lead to a fatal error 
%     (if the size of the empty times array is 0x0 rather than 0x1 or 1x0). 
%   - compute final bin correctly: In chronux 2.12, the final bin is (almost) always 0 if you
%     supply t, and 1/numTrials if you don't supply t. 
%   - note: new behavior with t not supplied uses final spike to set top
%     bin top; this means the final spike will not be counted
%   - note: does not remove undocumented final bin behavior when data is an array of
%     times and t is provided
% Usage: [dN,t]=binspikes(data,Fs,t)
% Inputs:
% data   (data as a structure array of spike times; or as a single
%        vector of spike times)
% Fs     (binning frequency)
% t      (the minimum and maximum times to be used to form the bins - [mint maxt]
%            - optional. Default use the spike times themselves to
%              determine the location of the bins. 
% Note: the times in data can be in any units. However, it is important
% that all units are chosen consistently. So, if spike times are in secs,
% Fs and t (if present) have to be in Hz and secs respectively. If spike
% times are in number of samples, Fs has to be 1, and t has to be in number
% of samples.
% Outputs:
% dN     (output binned spike counts as a matrix defined on bins starting with the
%         earliest spike across all channels and ending with the latest spike)
% t      (lower limit of each bin)
if nargin < 2; error('Need at least two input arguments'); end;
dt=1/Fs;
dtmp='';
if isstruct(data);
   C=length(data);
   fnames=fieldnames(data);
   if nargin <3 || isempty(t);
       mintime=zeros(1,C);
       maxtime=zeros(1,C);
       for ch=1:C
         eval(['dtmp=data(ch).' fnames{1} ';'])
         mintime(ch)=min(dtmp);
         maxtime(ch)=max(dtmp);
       end
       mintime=min(mintime);
       maxtime=max(maxtime)-dt; %modified to remove behavior that made last bin always 1/numTrials; sserene, 181022
   else
%        maxtimech=zeros(1,C);
%        for ch=1:C
%          eval(['dtmp=data(ch).' fnames{1} ';'])
% %          mintimech(ch)=min(dtmp);
%          maxtimech(ch)=max(dtmp);
%        end
       mintime=t(1);
       maxtime=t(end);
%        mintimech=min(mintimech);
%        maxtimech=max(maxtimech);
%        if maxtimech > max(t); t=[t maxtimech+dt]; end;
   end
   t=linspace(mintime,maxtime+dt,2+(maxtime-mintime)/dt); %modified to fix chronux bug that made last bin always zero; sserene, 180702
   dN = zeros(length(t)-1,C);
   for ch=1:C;
       eval(['dtmp=data(ch).' fnames{1} ';'])
       x=histc(dtmp,t);
       if ~isempty(x)            %modified to fix chronux bug that led to fatal errors when s.times was 0x0 
           dN(:,ch)= x(1:end-1); %modified to fix chronux bug that made last bin always zero; sserene, 180702
       end
   end
   t = t(1:end-1); %added to return correct bin bottoms with fix to chronux bug that made last bin always zero; sserene, 181022
else
   dtmp=data;
   if nargin < 3;
      mintime=min(dtmp);
      maxtime=max(dtmp)-dt; %modified to remove behavior that made last bin always 1/numTrials; sserene, 181022
   else
      mintime=t(1);
      maxtime=t(end);
   end
   t=linspace(mintime,maxtime+dt,2+(maxtime-mintime)/dt); %modified to fix chronux bug that made last bin always zero; sserene, 181022
%    if max(dtmp)>max(t); t=[t maxtime+dt]; end;
   x=histc(dtmp,t);
   if ~isempty(a)   %modified to fix chronux bug that led to fatal errors when data was 0x0 
     dN=x(1:end-1); %modified to fix chronux bug that made last bin always zero; sserene, 181022
   else
     dN = zeros(length(t)-1,1); %added to fix chronux bug that led to fatal errors when data was 0x0 
   end
   t = t(1:end-1); %added to return correct bin bottoms with fix to chronux bug that made last bin always zero; sserene, 181022
end



