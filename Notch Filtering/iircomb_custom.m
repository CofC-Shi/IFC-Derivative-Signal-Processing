function [num,den] = iircomb_custom(N,BW,varargin)

% THIS IS A MODIFIED VERSION OF THE STANDARD IIRCOMB() FUNCTION. I HAVE
% CHANGED THE VARIABLE maxLength ON LINE 59.

% Charlie Jindrich
% 05-30-24


%IIRCOMB IIR comb notching or peaking digital filter design.
%   [NUM,DEN] = IIRCOMB(N,BW) designs an Nth order comb notching digital 
%   filter with a -3 dB width of BW. N which must be a positive integer
%   is the number of notches in the range [0 2pi).  The notching filter 
%   transfer function is in the form
%
%                 1 - z^-N
%      H(z) = b * ---------
%                 1 - az^-N
%
%   The bandwidth BW is related to the Q-factor of a filter by BW = Wo/Q.
%
%   [NUM,DEN] = IIRCOMB(N,BW,Ab) designs a notching filter with a bandwidth
%   of BW at a level -Ab in decibels. If not specified, -Ab defaults to the 
%   -3 dB level (10*log10(1/2).  For peaking comb filters the default Ab is
%   3 dB or 10*log10(2).
%
%   [NUM,DEN] = IIRCOMB(...,TYPE) designs a comb filter with the specified
%   TYPE as either 'notch' or 'peak'.  If not specified it defaults to 'notch'.
%
%   The peaking filter transfer function is in the form
%
%                 1 + z^-N
%      H(z) = b * ---------
%                 1 - az^-N
%
%   EXAMPLE:
%      % Design a filter with a Q=35 to remove a 60 Hz periodic tone
%      % from system running at 600 Hz.
%      Fs = 600; Fo = 60;  Q = 35; BW = (Fo/(Fs/2))/Q;
%      [b,a] = iircomb(Fs/Fo,BW,'notch');  
%      fvtool(b,a);
% 
%   See also IIRNOTCH, IIRPEAK, FIRGR.

%   Author(s): P. Pacheco
%   Copyright 1999-2021 The MathWorks, Inc.

%   References:
%     [1] Sophocles J. Orfanidis, Introduction To Signal Processing
%         Prentice-Hall 1996.
%#codegen

narginchk(2,4);

% Validate input arguments.
[Ab,type] = combargchk(N,BW,varargin);

% Maximal allowed IIR filter order
maxLength = 2^10;
% Add range assertions, to preent codegen errors when DMA is disabled
assert(N <= maxLength);

if strncmp(type, 'notch',length(type))
    % Design a notching comb digital filter.
    [num,den] = notchingComb(N,BW,Ab);
else
    % Design a peaking comb digital filter.
    [num,den] = peakingComb(N,BW,Ab);
end

%------------------------------------------------------------------------
function [b,a] = notchingComb(N,BW,Ab)
% Design a comb digital filter.

% Inputs are normalized by pi.
BW = BW*pi;
D = N;

Gb   = 10^(-Ab/20);
beta = (sqrt(1-Gb.^2)/Gb)*tan(D*BW/4);
gain = 1/(1+beta);

ndelays = zeros(1,D-1);
b = gain*[1 ndelays -1];
a = [1 ndelays -(2*gain-1)];

%------------------------------------------------------------------------
function [b,a] = peakingComb(N,BW,Ab)
% Design a comb digital filter.

% Inputs are normalized by pi.
BW = BW*pi;
D = N;

Gb   = 10^(-Ab/20);
beta = (Gb/sqrt(1-Gb.^2))*tan(D*BW/4);
gain = 1/(1+beta);

ndelays = zeros(1,D-1);
b = (1-gain)*[1 ndelays -1];
a = [1 ndelays (2*gain-1)];

%------------------------------------------------------------------------
function [Ab,type] = combargchk(N,BW,opts)
% Checks the validity of the input arguments to IIRCOMB.

% Check the validity of order N.
validateattributes(N,{'single','double'},{'real','scalar','integer','>',0},'iircomb','N',1);

% Check the validity of BW.
validateattributes(BW,{'single','double'},{'real','scalar','>',0,'<',1},'iircomb','BW',2);

% Parse and validate optional input args.
[Ab,type] = parseoptions(opts);

%------------------------------------------------------------------------
function [Ab,type] = parseoptions(opts)
% Define default values.
AbDefault = abs(10*log10(.5)); % 3-dB width
typeDefault = 'notch';

% Parse the optional input arguments.  
switch length(opts)
    case 0
        % Define default values.
        Ab =  AbDefault; 
        type = typeDefault;
    case 1
        if ~dspIsChar(opts{1})
            Ab = checkAtten(opts{1});
            type = typeDefault;
        else
            type = checkCombType(char(opts{1})); % For comb filters.
            Ab = AbDefault;
        end
        
    case 2
        Ab = checkAtten(opts{1});
    	type = checkCombType(opts{2});    
end


%------------------------------------------------------------------------
function [Ab] = checkAtten(option)
% Determine if input argument is the attenuation scalar value.

validateattributes(option,{'single','double'},{'scalar'},'iircomb','Ab');
Ab = abs(option);

%------------------------------------------------------------------------
function [type] = checkCombType(option)
% Determine if input argument is the string 'notching' or 'peaking'.

type = validatestring(option,{'peak','notch'},'iircomb','type');

% [EOF]
