% IFC Signal Processing Comparison Script (Traditional vs Derivative Method)
clear; clc;

% --- Folder Settings for Multiple SNR Levels ---
noise_levels = {'5.5e-4', '3.5e-4', '2.2e-4', '1e-4', '5.5e-5', '3e-5', '2.2e-5'};
num_levels = length(noise_levels);
thresholds_trad = [1.5e-4, 1.5e-4, 1.5e-4, 1.5e-4, 5e-5,5e-5,5e-5];   % traditional method thresholds
thresholds_deriv = [2.5e-5, 2e-5, 2e-5, 2e-5, 1e-5, 1e-5, 1e-5];   % derivative method thresholds
snr_labels = [-5, 0, 5, 10, 15, 20, 25]; 

% --- Initialize Results ---
results_trad = struct('snr', [], 'comp_time', [], 'rms_err', [], 'event_acc', []);
results_deriv = results_trad;

% --- Process Each SNR Level ---
for n = 1:num_levels
    folder_path = "mat_data/" + noise_levels{n};
    clean_folder = folder_path + "_clean";
    file_list = dir(fullfile(folder_path, "*.mat"));

    fprintf("\nProcessing SNR level: %s\n", noise_levels{n});


    % --- Run Traditional Method ---
    [comp_time_trad, rms_err_trad, event_acc_trad] = run_traditional_method(folder_path, clean_folder, file_list,thresholds_trad(n));

    % --- Run Derivative Method ---
    [comp_time_deriv, rms_err_deriv, event_acc_deriv] = run_derivative_method(folder_path, clean_folder, file_list,thresholds_deriv(n));

    % --- Store Results ---
    results_trad(n).snr = noise_levels{n};
    results_trad(n).comp_time = mean(comp_time_trad);
    results_trad(n).rms_err = mean(rms_err_trad) * 1000;
    results_trad(n).event_acc = 100 - mean(event_acc_trad) * 100;
    

    results_deriv(n).snr = noise_levels{n};
    results_deriv(n).comp_time = mean(comp_time_deriv);
    results_deriv(n).rms_err = mean(rms_err_deriv) * 1000;
    results_deriv(n).event_acc = 100 - mean(event_acc_deriv) * 100;
    
end

% --- Plot Comparison ---

figure;
tiledlayout(1,3);

colors = [0.2 0.6 0.8; 0.6 0.6 0.6];  % [Trad color; Deriv color]

nexttile;
b = bar([ [results_trad.comp_time]', [results_deriv.comp_time]']);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
title('Computation Time'); ylabel('Time (s)');
% ylim([0 0.4]);
xticklabels(snr_labels); legend('Trad', 'Deriv'); grid on;

nexttile;
b = bar([ [results_trad.rms_err]', [results_deriv.rms_err]']);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
title('Amplitude RMS Error'); ylabel('Error (mV)');
% ylim([0 0.6]);
xticklabels(snr_labels); legend('Trad', 'Deriv'); grid on;

nexttile;
b = bar([ [results_trad.event_acc]', [results_deriv.event_acc]']);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
title('Event Detection Miss Rate'); ylabel('Missed Events (%)');
% ylim([0 1.0]);
xticklabels(snr_labels); legend('Trad', 'Deriv'); grid on;


%% --- Traditional Method Function ---
function [comp_time, avg_rms_err, event_accuracy] = run_traditional_method(folder_path, clean_folder, file_list, thresholds)
    num_files = length(file_list);
    comp_time = zeros(num_files, 1);
    event_accuracy = zeros(num_files, 1);
    avg_rms_err = zeros(num_files, 1);

    for i = 1:num_files
        tic;
        file_path = fullfile(folder_path, file_list(i).name);
        clean_path = fullfile(clean_folder, file_list(i).name);
        [signal, ~, Fs] = LoadData(file_path);
        sig = signal{2};  % channel 2
        sig = detrend(sig, 'linear');  % Remove baseline driftl

        gt_data = load(clean_path);
        x_clean = gt_data.dev18244.demods(2).sample.x;
        y_clean = gt_data.dev18244.demods(2).sample.y;
        sig_clean = sqrt(x_clean.^2 + y_clean.^2);
        sig_clean = detrend(sig_clean, 'linear');  % Remove baseline drift
        LM_gt = gt_data.LM_gt;
        RM_gt = gt_data.RM_gt;
        Map_gt = gt_data.Map;
        Mip_gt = gt_data.Mip;

        bpFilt = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',100,'HalfPowerFrequency2',1500,'SampleRate',Fs);
        filt_sig = filtfilt(bpFilt, sig);

        filt_sig = DataProcess(filt_sig, Fs);

        min_dist = round(0.5 * Fs / 1000);
        min_ht = thresholds;
        [pos_peaks, pos_locs] = findpeaks(filt_sig, 'MinPeakHeight', min_ht, 'MinPeakDistance', min_dist);
        [neg_peaks, neg_locs] = findpeaks(-filt_sig, 'MinPeakHeight', min_ht, 'MinPeakDistance', min_dist);
        neg_peaks = -neg_peaks;

        max_pp_samples = round(3.0 * Fs / 1000);
        j = 1;
        Time_Positive_ms = []; Time_Negative_ms = [];
        Amplitude_Positive = []; Amplitude_Negative = [];

        for k = 1:length(pos_locs)
            while j <= length(neg_locs) && (neg_locs(j) - pos_locs(k)) < 0
                j = j + 1;
            end
            if j > length(neg_locs), break; end
            dt = neg_locs(j) - pos_locs(k);
            if dt > 0 && dt <= max_pp_samples
                Time_Positive_ms(end+1) = (pos_locs(k)-1)/Fs*1000;
                Time_Negative_ms(end+1) = (neg_locs(j)-1)/Fs*1000;
                Amplitude_Positive(end+1) = pos_peaks(k);
                Amplitude_Negative(end+1) = neg_peaks(j);
                j = j + 1;
            end
        end

        %% Metrics
        % Event accuracy
        tolerance_samples = round(0.5 * Fs / 1000);  % e.g. ±0.5 ms tolerance
        TP = 0;
        
        for gt = LM_gt
            if any(abs(pos_locs - gt) <= tolerance_samples)
                TP = TP + 1;
            end
        end
        
        event_accuracy(i) = TP / length(LM_gt);

        % amplitude RMS Error
        match_pos = min(length(Map_gt), length(pos_peaks));
        match_neg = min(length(Mip_gt), length(neg_peaks));
        rms_pos = rms(pos_peaks(1:match_pos) - Map_gt(1:match_pos));
        rms_neg = rms(neg_peaks(1:match_neg) - Mip_gt(1:match_neg));
        avg_rms_err(i) = mean([rms_pos, rms_neg]);
        comp_time(i) = toc;

    end
end


%% --- Derivative Method Function ---
function [comp_time, avg_rms_err, event_accuracy] = run_derivative_method(folder_path, clean_folder, file_list, thresholds)
    num_files = length(file_list);
    comp_time = zeros(num_files, 1);
    event_accuracy = zeros(num_files, 1);
    avg_rms_err = zeros(num_files, 1);

    for i = 1:num_files
        tic;
        file_path = fullfile(folder_path, file_list(i).name);
        clean_path = fullfile(clean_folder, file_list(i).name);
        [signal, ~, Fs] = LoadData(file_path);
        sig = signal{2};  % channel 2
        sig = sig - mean(sig);  % DC offset removal

        gt_data = load(clean_path);
        x_clean = gt_data.dev18244.demods(2).sample.x;
        y_clean = gt_data.dev18244.demods(2).sample.y;
        sig_clean = sqrt(x_clean.^2 + y_clean.^2);
        sig_clean = sig_clean - mean(sig_clean);
        LM_gt = gt_data.LM_gt;
        RM_gt = gt_data.RM_gt;
        Map_gt = gt_data.Map;
        Mip_gt = gt_data.Mip;

        [reco, timing2, LM, RM, Map, Mip] = deriv_method(sig, thresholds, Fs);

        % Convert LM and RM to milliseconds
        Time_Positive_ms = (LM(:) - 1) / Fs * 1000; % column vector
        Time_Negative_ms = (RM(:) - 1) / Fs * 1000;
        Amplitude_Positive = Map(:);
        Amplitude_Negative = Mip(:);

        %% Accuracy Metrics
        % Event accuracy
        tolerance_samples = round(0.5 * Fs / 1000);  % e.g. ±0.5 ms tolerance
        TP = 0;
        
        for gt = LM_gt
            if any(abs(LM - gt) <= tolerance_samples)
                TP = TP + 1;
            end
        end
        
        event_accuracy(i) = TP / length(LM_gt);

        match_pos = min(length(Map_gt), length(Map));
        match_neg = min(length(Mip_gt), length(Mip));
        rms_pos = rms(Map(1:match_pos) - Map_gt(1:match_pos));
        rms_neg = rms(Mip(1:match_neg) - Mip_gt(1:match_neg));
        avg_rms_err(i) = mean([rms_pos, rms_neg]);
        comp_time(i) = toc;

    end
end
