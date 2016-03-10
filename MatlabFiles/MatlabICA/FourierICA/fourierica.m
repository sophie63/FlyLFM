%% FourierICA
% This simple implementation FourierICA is part of the supplementary material of the publication:
%
% A. Hyvärinen, P. Ramkumar, L. Parkkonen and R. Hari, _Independent component analysis of short-time Fourier transforms
% for spontaneous EEG/MEG analysis_, *NeuroImage* 49(1):257-271, 2010. 
%
% Should you use this code, we kindly request you to cite the aforementioned publication.
%
% <http://www.cs.helsinki.fi/group/neuroinf/code/fourierica/fourierica.m DOWNLOAD FourierICA from here> 
%
%

%% Overview
%
% FourierICA is an unsupervised learning method suitable for the analysis of rhythmic activity in EEG/MEG recordings.
% The method performs _independent component analysis_ (ICA) on short-time Fourier transforms of the data. As a result,
% more "interesting" sources with (amplitude modulated) oscillatory behaviour are uncovered and appropriately ranked.
% The method is also capable to reveal spatially-distributed sources appearing with different phases in different
% EEG/MEG channels by means of a complex-valued mixing matrix. For more details, please read the publication
%
% A. Hyvärinen, P. Ramkumar, L. Parkkonen and R. Hari, _Independent component analysis of short-time Fourier transforms
% for spontaneous EEG/MEG analysis_, *NeuroImage* 49(1):257-271, 2010,
%
% or go to the webpage <http://www.cs.helsinki.fi/u/ahyvarin/>
%
%

%% Input syntax formats
%
% # |[...]=fourierica(origdata,pcadim,samplfreq,minfreq,maxfreq)|
% # |[...]=fourierica(origdata,pcadim,samplfreq,minfreq,maxfreq,options{:})|
% # |[...]=fourierica(origdata,options{:})| 
% # |[...]=fourierica(origdata,options)|
%
% In (2) and (3) (but not in (4)), *options* is a cell list of the form {' _parname1_ ', _parvalue1_ , ...}. In (3),
% *options* must include the *mandatory parameters*.
%
% In (4), *options* is a structure whose fields have the same name as the parameters being specified.
%
% <html><h2><br></h2></html>

%% Output syntax
%
% # |*[S_FT,A_orig,W_orig]*=fourierica(origdata,varargin)| 
% # |*[S_FT,A_orig,W_orig,objective]*=fourierica(origdata,varargin)|
% # |*[S_FT,A_orig,W_orig,objective,convdata]*=fourierica(origdata,varargin)|
%

%% Input data
%
% Mandatory arguments are:
%
% # |*origdata*     :| EEG/MEG data of size channels x time points (rows x columns).
% # |*samplfreq*    :| Sampling frequency (in Hz).
% # |*pcadim*       :| Number of principal components to be extracted.
% # |*minfreq*      :| Minimum frequency to be analised (in Hz).
% # |*maxfreq*      :| Maximum frequency to be analised (in Hz).
%   
% Typical values would be: |pcadim=40,samplfreq=150, minfreq=7,maxfreq=30|.
%
% In addition, you can change any of the following optional parameters:
% 
% * |*complexmixing*    :| When |false| (default), mixing matrix is real-valued. Otherwise, complex-valued. 
% * |*windowlength_sec* :| Duration of short-time Fourier transform window (in seconds). By default, |windowlength_sec = 1|.
% * |*overlapfactor*    :| Overlap between consecutive windows. By default, |overlapfactor = 2|.
% * |*components*       :| Number of independent components. By default, |components = pcadim|.
% * |*hammingdata*      :| If |true|, data is multiplied by an appropriate Hamming window before short-time Fourier
% transformed. By default, |hammingdata = false|. More information on <matlab:doc('hamming') hamming>.
% * |*removeoutliers*   :| If |true| (default), remove outliers from data.
% * |*fcnoutliers*      :| Function handle of external code for removing outliers. By default, |fcnoutliers = 0|, meaning 
% that if |removeoutliers = 1|, then the internal removal takes place (see the code). To specify another function, write
% |fcnoutliers = @ _function-name_|.
% * |*maxiter*          :| Maximum number of iterations in estimating loop. By default, |maxiter = max(40*pcadim,2000)|.
% * |*conveps*          :| Convergence tolerance. By default, |conveps = 1e-7|.
% * |*zerotolerance*    :| Value of the smallest eigenvalue relative to the maximum eigenvalue that is considered
% different from zero. By default, |zerotolerance = conveps|.
% * |*showprogress*     :| if |true| (default), it shows progress information during the iteration. 
% * |*seed*             :| Sets the seed of the random number generator to a specific value. By default, |seed = -1|,
% which allows Matlab to choose the seed.
% * |*lambda*           :| Constant in objective function. 1 (default) or 0.1 seems to be the same.
%

%% Output data
%
% * S_FT: Source signals of the STFT's of data (3D tensor, channel x time x frequency). When mixing matrix is real, S_FT
% corresponds to the STFT's of the source signals in time domain, which can be obtained by computing W_orig*origdata.
% * A_orig: Spatial patterns in original space (not whitened)
% * W_orig: Spatial filters (pseudoinverse of A_orig)
% * objective: Values of the objective function for each independent component.
% * convdata: Structure with information of the convergence. Fields are: converged, convcriterion, numiter, pcadim and
% components.
%

%% Examples:
%
% * Example 1:
%   
%       [S,A,W]=fourierica(origdata,150,40,7,30);
%
% * Example 2:
%
%       options.samplfreq=150;
%       options.pcadim=40; 
%       options.minfreq=7; 
%       options.maxfreq=30; 
%       [S,A,W]=fourierica(origdata,options)
% 
% * Example 3:
%
%       [S,A,W]=fourierica(origdata,150,40,7,30,'complexmixing',true)
%
% * Example 4:
%
%       [S,A,W]=fourierica(origdata,'samplfreq',150,'pcadim',40,'minfreq',7,'maxfreq',30,'complexmixing',true)
%
% * Example 5:
%
%       options.samplfreq=150;
%       options.pcadim=40; 
%       options.minfreq=7; 
%       options.maxfreq=30; 
%       options.complexmixing=true;
%       [S,A,W]=fourierica(origdata,options) 
%
% The first two examples are equivalent. Likewise are the last three examples.
% Notice that the last three examples override the parameter _complexmixing_ , and thus complex mixing values
% are used instead.
%
% 

%% Version control
%
% *Original code* by Aapo Hyvärinen, University of Helsinki (8 Nov 2011).
% *Modified version* by Hugo Gabriel Eyherabide, University of Helsinki. (4 May 2012).
% 
% *Latest version* by Hugo Gabriel Eyherabide, University of Helsinki. (10 Sep 2013).
%
% Should you find a bug, please contact either Prof. Aapo Hyvärinen aapo.hyvarinen@helsinki.fi or 
% Hugo Gabriel Eyherabide hugo.eyherabide@helsinki.fi
%
%

%% License
%
% Copyright 2011 Aapo Hyvärinen
%
% Copyright 2012 Hugo Gabriel Eyherabide
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program.  If not, see
% <http://www.gnu.org/licenses/>.
%
% <html><h2><br></h2></html>




%% Explanation of the code
% <html><h2><br></h2></html>
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% HEADER
% </h2></li></ul></html>
%

function varargout = fourierica(origdata,varargin)

%%
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% INITIALIZATION
% </h2></li></ul></html>
%
% Input and output arguments are processed here. The code checks for available parameters and consistency. 
% Only parameters specified in parfixnames are mandatory. This parameters are all of type |numeric|.
% At the end, variables corresponding to the parameters are created in the function workspace.
%
fprintf('Launching Fourier-ICA\n');

% Check the number of output arguments
if nargout<3 || nargout>5, error('Wrong number of output arguments.'), end

% Determine number of input arguments 
numargin=size(varargin,2);

% Definition of mandatory parameters. They all are of type |numeric|. In the input, they must
% be introduced in the same order as in |parfixname|.
parfixnames={'pcadim','samplfreq','minfreq','maxfreq'};
% Additional mandatory parameters can be included by adding elements to the cell |parfixnames|. 

% Definition of optional parameters names (|paroptnames|), their types (|paropttypes|) and default values
% (|paroptdefval|). In the input, they can be introduced in any order using the input format shown in section "Input
% data"
paroptnames={'showprogress','seed','lambda','conveps','complexmixing','windowlength_sec','overlapfactor','components', ...
    'removeoutliers','maxiter','fcnoutliers','zerotolerance','hammingdata'};
paropttypes={'logical','numeric','numeric','numeric','logical','numeric','numeric','numeric',...
    'logical','numeric','function_handle','numeric','logical'};
paroptdefval={'true','-1','1','1e-7','0','1','2','pcadim','1','max(40*pcadim,2000)','[]','conveps','false'};
% Additional optional parameters can be included by adding elements to the cells |paroptnames|, |paropttypes| and
% |paroptdefval|. 


% Convert all input syntax formats into struct type 
numparfix=size(parfixnames,2);
numparopt=size(paroptnames,2);

% Check if there are enough input arguments.
if numargin>=numparfix
    % Check is mandatory parameters are of type |numeric|.
    if all(isnumeric([varargin{1:numparfix}]))
        for indpar=1:numparfix, eval(['param.' parfixnames{indpar} '=varargin{indpar};']), end
        argini=numparfix+2;
    else
        argini=2;
    end
    % Check if the number of optional parameters is correct.
    if rem(numargin-argini,2)==1, error('Optional parameters should always go by pairs'); end
    for indarg=argini:2:numargin, if ischar(varargin{indarg-1}), ...
                eval(['param.' varargin{indarg-1} '=varargin{indarg};']), end, end
elseif numargin==1
    param=varargin{1};
else
    error('Wrong number of input arguments.');
end

% Assign input arguments to function parameters (by creating variables in the function workspace).
numfields=0;

%
% Predefined messages for mandatory arguments
%
messages={'PCA dimension set to ';'Sampling frequency set to ';'Start of frequency band set to ' ...
    ;'End of frequency band set to '};

% Set mandatory parameters and create messages
for indpar=1:numparfix,
    if isfield(param,parfixnames{indpar}),
        if isnumeric(param.([parfixnames{indpar}]))
            eval([ parfixnames{indpar} '=param.(parfixnames{indpar});'])
            messages{indpar,2}=param.(parfixnames{indpar});
        else
            error([ parfixnames{indpar} ' should be numeric']);    
        end
        numfields=numfields+1;
    else
        error(['Field ' parfixnames{indpar} ' is mandatory.']);
    end
end

% Display messages about mandatory parameters settings
disp(messages);


% Cheak optional parameters type, set optional parameters and Set mandatory parameters and create messages
messages=cell(0,2);
for indpar=1:numparopt, 
    % Check if the optional parameter has been specified
    if isfield(param,paroptnames{indpar}), 
        numfields=numfields+1;
        % Check if the value of the optional parameter is of the correct type.
        if isa(param.(paroptnames{indpar}),paropttypes{indpar})
            % Assign value to optional parameter and creates message
            eval([ paroptnames{indpar} '=param.(paroptnames{indpar});'])
            messages{numfields-4,1}=['Further: ' paroptnames{indpar} ' set to '];
            messages{numfields-4,2}=param.(paroptnames{indpar});
        else
            error([ paroptnames{indpar} ' should be ' paropttypes{indpar}]);
        end
    else
        % Uses default value if it wasn't specified by the user.
        eval([ paroptnames{indpar} '=' paroptdefval{indpar} ';'])
    end
end

% Display messages for those optional parameters set according to the user specifications.
disp(messages);

% Warning about optional parameters not recognised.
if length(fieldnames(param))>numfields, ...
        warning('Some specified parameters were not recognized!'), end

clear param numargin varargin numfields paroptnames paroptdefval paropttypes parfixnames

%%
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% WINDOWING, SHORT-TIME FOURIER TRANSFORM AND FILTERING</h2></html>
% </h2></li></ul></html>
%
% The code first determines the correct window size and frequency band according to the parameters specified by the user
% and then computes the short-time Fourier transform (STFT). Filtering is implemented by an ideal band-pass filter
% (gain=1) between the frequencies specified by the user.
%

disp('Sampling windows and STFT');

% Determines number of channels (|channels|) and time points (|T|) in original data (|origdata|).
[channels,T]=size(origdata); 

%Compute number of time points in one window based on other parameters
windowsize=floor(windowlength_sec*samplfreq);
%Compute interval between consecutive windows
windowinterval=ceil(windowsize/overlapfactor);
%Compute number of windows
windows=floor((T-windowsize)/windowinterval+1);

%compute frequency indices (for the STFT)
startfftind=floor(minfreq*windowlength_sec+1); 
if startfftind<1, error('minfreq must be positive'); end
endfftind=floor(maxfreq*windowlength_sec+1); nyquistfreq=floor(windowsize/2);
if endfftind>nyquistfreq, error('maxfreq must be less than the Nyquist frequency'); end
fftsize=endfftind-startfftind+1;

% Initialization of tensor X, which is the main data matrix input to the code which follows.
X=zeros(fftsize,windows,channels);

% Define window initial limits
window=[1,windowsize];

% Construct Hamming window if necessary 
if hammingdata, hammingwindow=hamming(windowsize,'periodic'); end

% Short-time Fourier transform (window sampling + fft)
for j=1:windows
    % Extract data window
    datawindow=origdata(:,window(1):window(2));
    % Multiply by a Hamming window if necessary
    if hammingdata, datawindow=datawindow*diag(hammingwindow); end
    % Do FFT
    datawindow_ft=fft(datawindow');
    X(:,j,:)=datawindow_ft(startfftind:endfftind,:);
    % New data window interval
    window=window+windowinterval;
end

clear datawindow datawindow_ft window

%%
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% REMOVE OUTLIERS
% </h2></li></ul></html>
% 
% Outliers are defined as windows with large log-average power (_LAP_)
%
% $$LAP_{c,t}=log \sum_{f}{|X_{c,tf}|^2}$$
%
% where _c_, _t_ and _f_ are channels, window time-onsets and frequencies, respectively. The threshold is defined as
% |mean(LAP)+3 std(LAP)|. This process can be bypassed or replaced by specifying a function handle as an optional
% parameter.

if removeoutliers,
    if isempty(fcnoutliers)
        fprintf('Outlier removal:');
        lognorms=log(sum(abs(X.*conj(X)))); 
        outlierthreshold=mean(lognorms(:))+3*std(lognorms(:)); 
        outlierindices=find(lognorms>outlierthreshold);
        fprintf(' removed %u windows\n',length(outlierindices));
        X(:,outlierindices)=0; 
        clear lognorms outlierthreshold outlierindices
    else
        X=fcnoutliers(X); 
    end
end


%%
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% PCA, DIMENSION REDUCTION, WHITENING
% </h2></li></ul></html>
%
% Do principal component analysis, reduce dimension of the data and whitening with respect to channels. This is
% compulsory, as typical in BSS algorithms.
% 

disp('Spatially whitening data with PCA dimension reduction');
% Concatenate STFT for consecutive windows in each channel
Xmat_c=reshape(X,[fftsize*windows,channels]);
% Store sample size after this transformation
N=fftsize*windows;
clear X

% Substract mean value from channels
Xmat_c_mv=mean(Xmat_c);
Xmat_c=Xmat_c-Xmat_c_mv(ones(N,1),:);
clear Xmat_c_mv

%Do PCA (on matrix Xmat_c)
covmat=cov(Xmat_c);
if ~complexmixing, covmat=real(covmat); else covmat=conj(covmat); end
[Ec, Dc] = eig(covmat);
[d,order] = sort(diag(Dc),'descend');
clear covmat Dc

% Checks for negative eigenvalues
if any(d(1:pcadim)<0), warning('Negative eigenvalues! Reducing PCA and ICA dimension...'), end
% Check for eigenvalues near zero (relative to the maximum eigenvalue)
zeroeigval=sum((d(1:pcadim)/d(1))<zerotolerance);
% Adjust dimensions if necessary (because zero eigenvalues were found)
pcadim=pcadim-zeroeigval;
if pcadim<components, components=pcadim; end
if zeroeigval, fprintf('PCA dimension is %d and ICA dimension is %d\n',pcadim,components); end


% Construct whitening and dewhitening matrices
dsqrt = sqrt(d(1:pcadim));
dsqrtinv = 1./dsqrt;
Ec = Ec(:,order(1:pcadim));
whiteningmatrix=diag(dsqrtinv)*Ec';
dewhiteningmatrix=Ec*diag(dsqrt);
clear d order zeroeigval dsqrt dsqrtinv Ec


% Reduce dimensions and whiten data. |Zmat_c| is the main input for the iterative algorithm
Zmat_c=whiteningmatrix*transpose(Xmat_c);
Zmat_c_tr=Zmat_c'; % Also used in the fixed-point iteration.
clear Xmat_c



%%
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% COMPLEX-VALUED FAST-ICA ESTIMATION
% </h2></li></ul></html>
%

fprintf('Launching complex-valued FastICA:                              ');

% Sets the random number generator
switch seed
    case -1, numrandn=@randn;
    otherwise, stream=RandStream('mrg32k3a','seed',19); numrandn=@(x,y)randn(stream,x,y);
end


% Initial point, make it imaginary and unitary
switch complexmixing, 
    case true, W_old=complex(numrandn(components,pcadim),numrandn(components,pcadim)); 
    case false, W_old=numrandn(components,pcadim); 
end
W_old=sqrtm(inv(W_old*W_old'))*W_old;


% Iteration starts here
  
tic
for iter=1:maxiter
    %Compute outputs, note lack of conjugate
    Y=W_old*Zmat_c;
    
    %Computing nonlinearities
    Y2=abs(Y.*conj(Y));
    gY2 = 1./(lambda + Y2);
    dmv=lambda*sum(gY2.^2,2);
    
    %Fixed-point iteration
    W_new=(Y.*gY2)*Zmat_c_tr-diag(dmv)*W_old;
    
    %In case we want to restrict W to be real-valued, do it here:
    switch complexmixing, case true, W=W_new; case false, W=real(W_new); end

    %Make unitary
    W=sqrtm(inv(W*W'))*W;
    
    %check if converged
    convcriterion=1-sum(abs(sum(W.*conj(W_old),2)))/components;
    if convcriterion<conveps, break, end
    
    %Show progress to user
    if showprogress, 
        erasechar={'\b'};
        infostr=[num2str(iter,'%4d') '    ' num2str(convcriterion,'%1.4e')];
        fprintf([erasechar{ones(length(infostr),1)} infostr]);
   end
    
    %store old value
    W_old=W;
end
toc

clear Y Y2 gY2 dmv W_old W_new Zmat_c_tr

%Compute mixing matrix (in whitened space)
A=W';
%Compute source signal estimates
S=W*Zmat_c;

%Tell if convergence problems
if convcriterion>conveps,
    fprintf(2,'\nFailed to converge, results may be wrong!\n')
else
    fprintf('\nConverged.\n')
end


%%
%
% <html><ul style="margin-left:-1.5em;"><li><h2 style="color:#000077;">
% SORT COMPONENTS AND TRANSFORMS TO ORIGINAL SPACE
% </h2></li></ul></html>
%

fprintf('Sorting components and reformatting results.\n')

%compute objective function for each component
objective=-mean(log(lambda+abs(S.*conj(S))),2);

%sort components using the objective
[objective,componentorder]=sort(objective,'descend');
W=W(componentorder,:);
A=A(:,componentorder);
S=S(componentorder,:);

%Compute mixing and demixing matrix in original channel space
%Spatial filters
W_orig=W*whiteningmatrix;
%Spatial patterns
A_orig=dewhiteningmatrix*A;
%Independent components. If real mixing matrix, these are STFT's of the
%sources.
S_FT=permute(reshape(S,[components,fftsize,windows]),[1,3,2]);
%Output convergence data

switch nargout
    case 3, varargout={S_FT,A_orig,W_orig};
    case 4, varargout={S_FT,A_orig,W_orig,objective};
    case 5, convdata=struct('convcriterion',convcriterion,'converged',convcriterion<conveps, ...
            'numiter',numiter,'pcadim',pcadim,'components',components); ...
            varargout={S_FT,A_orig,W_orig,objective,convdata};
end

fprintf('Done.\n')


 
