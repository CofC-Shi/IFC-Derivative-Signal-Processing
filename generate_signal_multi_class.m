function [r_noisy, r_clean, LM_gt, RM_gt, Map, Mip, true_labels] = generate_signal_multi_class(Fs, bead_diameters, bead_labels)
% generate_signal_multi_class - Generate synthetic IFC signals for multiple classes
% Outputs include lock-in demodulated signal and ground truth event positions for training

% Inputs:
%   Fs             - Sampling frequency (Hz)
%   bead_diameters - Array of bead diameters (m)
%   bead_labels    - Class labels for each bead diameter
%
% Outputs:
%   r            - Demodulated lock-in magnitude signal
%   LM_gt, RM_gt - Ground truth left/right midpoints of events
%   Map, Mip     - Event peak amplitudes
%   true_labels  - Ground truth class labels for each event

    %% Parameters
    duration = 30;                  % seconds
    carrier_freq = 1.1e6;           % Hz
    num_points = round(Fs * duration);
    offset_value = 1e-2;
    pink_noise_level = 5e-6;   % default 5e-6;
    white_noise_level = 5e-7;  % default 5e-7
    periodic_noise_level = 2.2e-5;  % default 2.2e-5
    lp_cutoff =500;               % low-pass cutoff after demodulation
    event_rate = 5;                % events/sec

    noise_freqs = [3.5, 6.6, 10, 20, 26.7, 30, 40, 50, 60, 70, 80, ...
        90, 100, 110, 120, 130, 140, 150, 160, 170, 180];

    noise_amps = [0.0658684, 0.2805874, 0.4893336, 0.0592756, ...
        0.1324801, 0.1, 0.06363636, 0.06818182, 1.0, 0.04045455, ...
        0.03772727, 0.04363636, 0.04272727, 0.03454545, ...
        0.07727273, 0.02590909, 0.03636364, 0.02454545, ...
        0.02090909, 0.01636364, 0.05909091] * periodic_noise_level;

    t = (0:num_points-1) / Fs;

    %% Initialize signals
    signal = offset_value * ones(1, num_points);
    LM_gt = []; RM_gt = []; true_labels = [];
    Map = []; Mip = [];

    %% Generate Poisson-distributed events
    expected_events = round(event_rate * duration);
    inter_arrival = exprnd(1/event_rate, 1, expected_events);
    event_times = cumsum(inter_arrival);
    event_times = event_times(event_times < duration);
    event_locs = round(event_times * Fs);

    for i = 1:length(event_locs)
        idx = randi(length(bead_diameters));
        d_mean = bead_diameters(idx);
        d_sigma = 0.05 * d_mean;  % 5% size variation
        d = abs(normrnd(d_mean, d_sigma));  % ensure positive size

        label = bead_labels(idx);

        amp = 2 * abs((d / 1.1e-4)^3);  % scaling factor
        transit =  (-293.4 * d + 0.004)/20;
        %disp(transit);
        width = 0.3 * transit;

        tc = event_locs(i);
        if tc <= 1 || tc >= num_points
            continue;
        end

        ex1 =  amp * exp(-((t - t(tc)) + transit/2).^2 / (2 * width^2));
        ex2 = -amp * exp(-((t - t(tc)) - transit/2).^2 / (2 * width^2));
        signal = signal + (ex1 + ex2);

        % Ground truth bookkeeping
        LM_gt = [LM_gt, tc - (transit/2 * Fs)];
        RM_gt = [RM_gt, tc + (transit/2 * Fs)];
        Map = [Map, amp];
        Mip = [Mip, -amp];
        true_labels = [true_labels, label];
    end

    %% Lock-in Modulation and Demodulation
    carrier = cos(2 * pi * carrier_freq * t);
    modulated_signal = signal .* carrier;

    ref_cos = cos(2*pi*carrier_freq*t);
    ref_sin = sin(2*pi*carrier_freq*t);
    x_raw = modulated_signal .* ref_cos;
    y_raw = modulated_signal .* ref_sin;

    %–– Design a 4‑th‑order Butterworth LPF
    Wn   = lp_cutoff/(Fs/2);          % normalized cutoff (0–1, Nyquist‑scaled)
    [b,a] = butter(4, Wn, 'low');     % 4th‑order Butterworth coefficients
    
    %–– Apply with filtfilt for zero‑phase (forward + reverse)
    x_demod = filtfilt(b, a, x_raw);
    y_demod = filtfilt(b, a, y_raw);

    r_clean = sqrt(x_demod.^2 + y_demod.^2);

    %% Add noise after demodulation
    periodic_noise = zeros(1, num_points);
    for j = 1:length(noise_freqs)
        phi = 2*pi*rand;
        periodic_noise = periodic_noise + ...
            (noise_amps(j)/2) * sin(2*pi*noise_freqs(j)*t + phi);
    end

    pink = generate_pink_noise(num_points, Fs) * pink_noise_level;
    x_demod = x_demod + pink + periodic_noise;
    y_demod = y_demod + pink + periodic_noise;

    % --- Generate white noise ---
    white_noise = randn(1, num_points);                    % Standard white noise
    white_noise = white_noise / max(abs(white_noise));     % Normalize
    white_noise = white_noise * white_noise_level;         % Scale to desired amplitude

    % --- Add to demodulated signals ---
    x_demod = x_demod + white_noise;
    y_demod = y_demod + white_noise;
    
    % --- Final noisy signal ---
    r_noisy = sqrt(x_demod.^2 + y_demod.^2);

    % === Debugging Visualization ===
    debug_plot = false;
    
    if debug_plot
        figure;
        t_ms = t * 1000;
        subplot(2,1,1); plot(t_ms, r_clean); title('Original Bipolar Signal'); ylabel('V'); grid on;
        xlim([10 1000]);
        subplot(2,1,2); plot(t_ms, r_noisy); title('Recovered Magnitude from X/Y (Lock-in)'); xlabel('Time (ms)'); ylabel('|Z|'); grid on;
        xlim([10 1000]);
    end
end
   

%% Pink noise helper
function pink_noise = generate_pink_noise(N, Fs)
    white = randn(1, N);
    F = fft(white);
    f = (0:N/2) * Fs / N; f(1) = 1;
    filter = 1 ./ sqrt(f);
    F_filtered = F(1:N/2+1) .* filter;
    F_full = [F_filtered, conj(F_filtered(end-1:-1:2))];
    pink_noise = real(ifft(F_full));
    pink_noise = pink_noise / max(abs(pink_noise));
end
