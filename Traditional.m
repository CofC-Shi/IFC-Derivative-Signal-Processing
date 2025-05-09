%% IFC Signal Processing Comparison Script (Traditional vs Derivative Method)
clear; clc;

% --- Folder Paths ---
folder_path = "mat_data/fast_1.0";
clean_folder = folder_path + "_clean";
file_list = dir(fullfile(folder_path, "*.mat"));
num_files = length(file_list);

fprintf("Found %d files for processing in '%s'\n", num_files, folder_path);


%% --- Traditional Method Function ---

% --- Initialize Metrics ---
comp_time = zeros(num_files, 1);
event_accuracy = zeros(num_files, 1);
transit_accuracy = zeros(num_files, 1);
avg_rms_err = zeros(num_files, 1);
TH = 2.5e-5;


for i = 1:1%length(file_list) % Change 1 to full batch when needed
    tic; % Start timer

    disp(['Processing file ', num2str(i), ' of ', num2str(length(file_list)), ': ', file_list(i).name]);
    % --- Load Signal ---
    file_path = fullfile(folder_path, file_list(i).name);
    [signal, labels, Fs] = LoadData(file_path); 
    sig = signal{2};  % channel 2
    sig = detrend(sig, 'linear');  % Remove baseline drift

    % --- Load Clean Signal ---
    % Extract the name part (e.g., '1e-5') and build the clean path
    clean_path = fullfile(clean_folder, file_list(i).name);

    disp(['Clean file path: ', clean_path]);  % <-- Print the clean path
    clean_data = LoadData(clean_path);
    sig_clean = clean_data{2};
    sig_clean = detrend(sig_clean, 'linear');

    % --- Bandpass Filter (0.1–1.5 kHz typical range) ---
    bpFilt = designfilt('bandpassiir', 'FilterOrder', 4, ...
             'HalfPowerFrequency1', 100, 'HalfPowerFrequency2', 1500, ...
             'SampleRate', Fs);
    filt_sig = filtfilt(bpFilt, sig);

    % % Notch filter and detrend
    % filt_sig = DataProcess(filt_sig, Fs);

    % --- Detect Peaks Directly on Filtered Signal ---
    min_dist = round(0.5 * Fs / 1000);
    min_ht = TH;

    [pos_peaks, pos_locs] = findpeaks(filt_sig, 'MinPeakHeight', min_ht, 'MinPeakDistance', min_dist);
    [neg_peaks, neg_locs] = findpeaks(-filt_sig, 'MinPeakHeight', min_ht, 'MinPeakDistance', min_dist);
    neg_peaks = -neg_peaks;

    % --- Align Pairs with Max Time Constraint ---
    max_pp_interval_ms = 5.0;
    max_pp_samples = round(max_pp_interval_ms * Fs / 1000);

    Time_Positive_ms = [];
    Time_Negative_ms = [];
    Amplitude_Positive = [];
    Amplitude_Negative = [];

    j = 1;
    for k = 1:length(pos_locs)
        while j <= length(neg_locs) && (neg_locs(j) - pos_locs(k)) < 0
            j = j + 1; % Skip neg before pos
        end
        if j > length(neg_locs), break; end

        dt = neg_locs(j) - pos_locs(k);
        if dt > 0 && dt <= max_pp_samples
            % Valid pair
            Time_Positive_ms(end+1) = (pos_locs(k)-1)/Fs*1000;
            Time_Negative_ms(end+1) = (neg_locs(j)-1)/Fs*1000;
            Amplitude_Positive(end+1) = pos_peaks(k);
            Amplitude_Negative(end+1) = neg_peaks(j);
            j = j + 1;  % Move to next negative peak
        end
    end


    metric = true;
    if metric
        % --- Transit Time Accuracy (Positive to Negative Peak)
        % Make sure peaks are paired
        num_pairs = min(length(Time_Positive_ms), length(Time_Negative_ms));
        transit_pred = Time_Negative_ms(1:num_pairs) - Time_Positive_ms(1:num_pairs);

        events = detect_ev(sig_clean, 0);
        gt_transit = [];
        for k = 1:2:length(events)-1
            seg = sig_clean(events(k):events(k+1));
            [pks_pos, locs_pos] = findpeaks(seg);
            [pks_neg, locs_neg] = findpeaks(-seg);
            if ~isempty(locs_pos) && ~isempty(locs_neg)
                t_pos = (events(k) + locs_pos(1) - 1) / Fs * 1000; % ms
                t_neg = (events(k) + locs_neg(1) - 1) / Fs * 1000; % ms
                gt_transit(end+1) = t_neg - t_pos;
            end
        end

        match_len = min(length(gt_transit), length(transit_pred));
        if match_len >= 1
            transit_error = abs(transit_pred(1:match_len) - gt_transit(1:match_len));
            transit_accuracy(i) = mean(transit_error) / mean(gt_transit(1:match_len));
        else
            transit_accuracy(i) = NaN;
        end

        % --- Accuracy based on clean signal event windows ---
        true_events = length(events) / 2;
        detected_events = length(Time_Positive_ms);
        event_accuracy(i) = detected_events / true_events;

        % --- RMS Error (positive and negative peaks) ---
        [clean_pos, ~] = findpeaks(sig_clean);
        [clean_neg, ~] = findpeaks(-sig_clean); clean_neg = -clean_neg;

        match_pos = min(length(clean_pos), length(pos_peaks));
        match_neg = min(length(clean_neg), length(neg_peaks));
        rms_pos = rms(pos_peaks(1:match_pos) - clean_pos(1:match_pos));
        rms_neg = rms(neg_peaks(1:match_neg) - clean_neg(1:match_neg));
        avg_rms_err(i) = mean([rms_pos, rms_neg]);

        comp_time(i) = toc;
    end

    % --- Plot Results ---
    plotting = true;
    if plotting
        t = (0:length(sig)-1) / Fs * 1000;
        figure;
        tiledlayout(2,1);

        % Original Signal
        nexttile;
        plot(t, sig, 'Color', [0.6 0.6 0.6]);hold on;
        plot(t, sig_clean(1:length(sig)), 'g');  % Green = clean signal
        xlabel('Time (ms)');
        ylabel('Amplitude (V)');
        title('Original Signal');
        legend('Noisy Signal', 'Clean Signal');
        grid on;

        % Filtered Signal
        nexttile;
        plot(t, filt_sig, 'b'); hold on;

        % Plot detected peaks
        scatter(Time_Positive_ms, Amplitude_Positive, 'r', 'filled');
        scatter(Time_Negative_ms, Amplitude_Negative, 'b', 'filled');

        % Connect positive to negative peaks
        num_pairs = min(length(Time_Positive_ms), length(Time_Negative_ms));
        for k = 1:num_pairs
            plot([Time_Positive_ms(k), Time_Negative_ms(k)], ...
                 [Amplitude_Positive(k), Amplitude_Negative(k)], ...
                 '--k', 'LineWidth', 1.2);
        end
        xlabel('Time (ms)');
        ylabel('Amplitude (V)');
        title('Filtered Signal');
        legend('Filtered Signal', 'Positive Peak', 'Negative Peak', 'Connection');
        grid on;
    end
end

% --- Summary ---
fprintf("Average Computation Time: %.4f s\n", mean(comp_time));
fprintf("Average RMS Error: %.2e V\n", mean(avg_rms_err));
fprintf("Average Event Accuracy: %.2f %%\n", mean(event_accuracy) * 100);
fprintf("Average Transit Time Accuracy: %.2f %%\n", mean(transit_accuracy) * 100);
