% add boundary checking by LS on 11/07/2024

function [LL, LM, RM, RR] = LaR(px, sig)
% LaR - Identifies left and right boundaries around each peak in a signal
% using zero-crossing points.
%
% Inputs:
%   px  - Array of peak indices in the signal
%   sig - Signal data array
%
% Outputs:
%   LL - Left boundary zero-crossing index for each peak
%   LM - Left midpoint zero-crossing index for each peak
%   RM - Right midpoint zero-crossing index for each peak
%   RR - Right boundary zero-crossing index for each peak

% Find zero-crossing indices of the signal
zero_crossings = zeroz(sig);
num_zero_crossings = length(zero_crossings);

% Initialize output arrays
num_peaks = length(px);
LL = zeros(1, num_peaks);
LM = zeros(1, num_peaks);
RM = zeros(1, num_peaks);
RR = zeros(1, num_peaks);

% Loop through each peak
for k = 1:num_peaks
    % Find the first zero-crossing index greater than the current peak
    ix = find(zero_crossings > px(k), 1);
    
    % Ensure we have enough zero-crossings around the peak
    if isempty(ix) || ix < 3 || ix > num_zero_crossings
        continue; % Skip if index is out of bounds or there aren’t enough surrounding zero-crossings
    end
    
    % Assign left and right midpoints around the peak
    LM(k) = zero_crossings(ix - 1);
    RM(k) = zero_crossings(ix);

    % Initialize indices and counters for boundary search
    zx = ix - 2; % Start left search before LM
    xz = ix + 1; % Start right search after RM
    max_iterations = 50;
    
    % Find left boundary (LL), ensuring sufficient spacing from LM
    count = 0;
    while zx > 0 && (LM(k) - zero_crossings(zx)) < round(0.5 * (RM(k) - LM(k))) && count < max_iterations
        zx = zx - 1;
        count = count + 1;
    end
    LL(k) = zero_crossings(max(zx, 1)); % Assign valid LL or default to closest valid index
    
    % Find right boundary (RR), ensuring sufficient spacing from RM
    count = 0;
    while xz <= num_zero_crossings && (zero_crossings(xz) - RM(k)) < round(0.5 * (RM(k) - LM(k))) && count < max_iterations
        xz = xz + 1;
        count = count + 1;
    end
    RR(k) = zero_crossings(min(xz, num_zero_crossings)); % Assign valid RR or default to closest valid index
end

