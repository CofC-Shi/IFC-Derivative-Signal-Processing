% Derivative Method for Signal Processing
% This function processes the input signal using a derivative-based method to identify events,
% reconstructs the signal, and measures the execution time.

function [reco, timing2, LM, RM, Map, Mip] = deriv_method(sig, TH, Fs)

% Start timing the execution
t_start = clock;

% % % Calculate the derivative of the signal using a custom BDC function
deriv = BDC(sig);

% % % Define sampling parameters
% num_samples = length(sig);
% time_vector = (0:num_samples-1) / Fs; % Time vector in seconds
% 
% % Plot the Derivative of the Signal
% plot(time_vector, deriv * 1000); % Convert to mV/s for visualization
% title('Derivative of Signal');
% xlabel('Time (s)');
% ylabel('Derivative (mV/s)');

% Detect Peaks in the Derivative Signal
[px, py] = NP(TH, deriv);

% Find Zero Crossings for Labeling Events
[LL, LM, RM, RR] = LaR(px, deriv);

% Calculate Peak Values from Detected Events
[Map, Mip] = calcpeak(LL, LM, RM, RR, sig);

% --- Remove Duplicate Events (based on LM and RM pairs) ---
pairs = [LM(:), RM(:)];
[~, unique_idx] = unique(pairs, 'rows', 'stable');

LM = LM(unique_idx);
RM = RM(unique_idx);
Map = Map(unique_idx);
Mip = Mip(unique_idx);

% Reconstruct the Signal using a Bi-Gaussian Model
reco = bigaus_recon(LM, RM, Map, Mip, deriv);


% Measure the Execution Time in Milliseconds
timing2 = round(1000 * etime(clock, t_start));

end


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
    if isempty(ix) || ix < 2 || ix > num_zero_crossings
        continue; % Skip if index is out of bounds or there aren’t enough surrounding zero-crossings
    end

    % Assign left and right midpoints around the peak
    LM(k) = zero_crossings(ix - 1);
    RM(k) = zero_crossings(ix);

    % Initialize indices and counters for boundary search
    zx = ix - 2; % Start left search before LM
    xz = ix + 1; % Start right search after RM
    max_iterations = 50;

    half_width = round(RM(k) - LM(k));
    min_dist = round(0.5 * half_width);  % Minimum spacing requirement

    % ----- Left boundary (LL) -----
    count = 0;
    LL_found = false;
    while zx > 0 && (LM(k) - zero_crossings(zx)) >= min_dist && (LM(k) - zero_crossings(zx)) <= 1.5 * half_width && count < max_iterations
        LL(k) = zero_crossings(zx);
        LL_found = true;
        break;
    end
    if ~LL_found
        LL(k) = LM(k) - half_width;  % fallback with guaranteed minimum spacing
    end

    % ----- Right boundary (RR) -----
    count = 0;
    RR_found = false;
    while xz <= num_zero_crossings && (zero_crossings(xz) - RM(k)) >= min_dist && (zero_crossings(xz) - RM(k)) <= 1.5 * half_width && count < max_iterations
        RR(k) = zero_crossings(xz);
        RR_found = true;
        break;
    end
    if ~RR_found
        RR(k) = RM(k) + half_width;  % fallback
    end
end

% --- Fix overlapping RR and LL by setting both to midpoint of RM and LM ---
for k = 1:num_peaks - 1
    if RR(k) >= LL(k+1)
        % Compute midpoint between RM(k) and LM(k+1)
        split_point = round((RM(k) + LM(k+1)) / 2);

        % Update both boundary points to this split point
        RR(k) = split_point;
        LL(k+1) = split_point;
    end
end


% ===== Debug Plot =====
debug_plot = false;  % Set to true to show debug figure

if debug_plot
    figure;
    t = 1:length(sig);
    plot(t, sig, 'k', 'LineWidth', 1.2); hold on;

    % Only plot valid (non-zero) points
    valid_LL = LL(LL > 0);
    valid_LM = LM(LM > 0);
    valid_RM = RM(RM > 0);
    valid_RR = RR(RR > 0);

    scatter(valid_LL, sig(round(valid_LL)), 'go', 'filled');    % Left boundary (Green)
    scatter(valid_LM, sig(round(valid_LM)), 'bo', 'filled');     % Left midpoint (Blue)
    scatter(valid_RM, sig(round(valid_RM)), 'ro', 'filled');     % Right midpoint (Red)
    scatter(valid_RR, sig(round(valid_RR)), 'mo', 'filled');     % Right boundary (Magenta)

    legend('Signal', 'LL (Green)', 'LM (Blue)', 'RM (Red)', 'RR (Magenta)', 'Location', 'best');
    xlim([0 9000]);
    xlabel('Sample Index');
    ylabel('Amplitude');
    title('Debug plot of LL, LM, RM, RR Points');
    grid on;
    hold off;
end


end


% find the peaks and return timestamp (peak indices) and peak amplitude

function [x, y] = NP(TH, sig)

    % Take derivative
    sigD = BDC(sig);

    % Find zero-crossings where derivative goes from negative to positive (minima)
    N2P = find(sigD(1:end-1) < 0 & sigD(2:end) > 0);

    % Parameters
    min_peak_dist = 5;  % Adjust as needed based on expected event width

    x = [];
    y = [];

    for k = 1:length(N2P)
        idx = N2P(k);
        val = sig(idx);

        % Threshold condition
        if val <= -TH
            if isempty(x)
                x(end+1) = idx;
                y(end+1) = val;
            else
                % Check distance from previous peak
                if (idx - x(end)) >= min_peak_dist
                    x(end+1) = idx;
                    y(end+1) = val;
                elseif val < y(end)  % Prefer deeper minimum
                    x(end) = idx;
                    y(end) = val;
                end
            end
        end
    end
end


function [Map, Mip] = calcpeak(LL, LM, RM, RR, sig)
    deriv = BDC(sig);

    Map = zeros(1, length(LL));
    Mip = zeros(1, length(LL));

    for k = 1:length(LL)
        if LM(k) == 0 || RM(k) == 0
            continue
        end

        % Round indices properly
        ll_idx = round(LL(k));
        lm_idx = round(LM(k));
        rm_idx = round(RM(k));

        % Avoid invalid indexing
        ll_idx = max(1, min(length(deriv), ll_idx));
        lm_idx = max(1, min(length(deriv), lm_idx));
        rm_idx = max(1, min(length(deriv), rm_idx));
        
        Map(k) = sum(deriv(ll_idx:lm_idx)) + sig(ll_idx);
        Mip(k) = sum(deriv(ll_idx:rm_idx)) + sig(ll_idx);

    end
end


function x = zeroz(sig)

%peaks
P2N = find(sig(1:end-1) > 0 & sig(2:end) < 0);
%mins
N2P = find(sig(1:end-1) < 0 & sig(2:end) > 0);
x = sort([P2N, N2P]);
end