function [processed_data,timing2] = DataProcess(raw_data, sample_freq)

t2 = clock;

% DESCRIPTION
%
%   Filters and detrends raw data. 60 Hz and its primary harmonics (120 and
%   180) are filtered using notch filters. A comb filter is used to remove
%   10 Hz plus its harmonics. The signal is detrended at the end to remove
%   any signal drift and any non-linearities introduced by the filtering.
%
%   Author - Charlie Jindrich
%   Date - 05/09/24

% INPUTS
%
%   raw_data - vector containing raw impedance data.
%
%   sample_freq - Sampling frequency.

% OUTPUTS
%  
%   processed_data - vector of processed data.

%STEP 1
L = length(raw_data);

% Breakpoint Vector
res = 200;
num_steps = round(L/res);
bp = (1:num_steps)*res;


%timing = round(1000 * etime(clock, t2));
%disp(['Step 1 time: ',num2str(timing),' miliseconds']);


%STEP 2

% Notch filter Parameters
notch_freq_1 = 60; % Frequency to remove (Hz)
bandwidth_1 = 3; % Bandwidth around the notch frequency
[bn1, an1] = iirnotch(notch_freq_1/(sample_freq/2), bandwidth_1/(sample_freq/2));

notch_freq_2 = 120; % Frequency to remove (Hz)
bandwidth_2 = 3; % Bandwidth around the notch frequency
[bn2, an2] = iirnotch(notch_freq_2/(sample_freq/2), bandwidth_2/(sample_freq/2));

notch_freq_3 = 180; % Frequency to remove (Hz)
bandwidth_3 = 3; % Bandwidth around the notch frequency
[bn3, an3] = iirnotch(notch_freq_3/(sample_freq/2), bandwidth_3/(sample_freq/2));

% Define Comb Filter
q = 20;     % Scales bandwidth
fo = 10;    % Freq to remove
bw = (fo/(sample_freq/2))/q;


%timing = round(1000 * etime(clock, t2));
%disp(['Step 2 time: ',num2str(timing),' miliseconds']);

%Step 3

% iircomb_custom is a copy of the standard iircomb() function except that I
% changed the maximum allowable filter order.
[bn4,an4] = iircomb_custom(round(sample_freq/fo),bw); % Note type flag 'notch'

r = raw_data;

%-------------------------------- Filter

r = filtfilt(bn1, an1, r);
r = filtfilt(bn2, an2, r);
r = filtfilt(bn3, an3, r);
r = filtfilt(bn4, an4, r);

%-------------------------------- Detrend
%    timing = round(1000 * etime(clock, t2));
%disp(['Step 3 time: ',num2str(timing),' miliseconds']);

%Step 4
%r = detrend(r, 1, bp);
r = Linear_Detrend(r,100);

processed_data = r;

%timing = round(1000 * etime(clock, t2));
%disp(['Detrend time: ',num2str(timing),' miliseconds']);

timing2 = round(1000 * etime(clock, t2));
%disp(['Notch Total time: ',num2str(timing2),' miliseconds']);

end
