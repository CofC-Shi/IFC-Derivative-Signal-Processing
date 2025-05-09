
clear all; close all; clc;

Fs = 7196;         
TH = 5e-6;        
bead_diameters = [6e-6, 5e-6, 4e-6];
bead_labels = [1, 2, 3];

stream_count = 30;
all_features = [];
all_labels = [];
fprintf('\n Retraining classifier with %d simulated streams...\n', stream_count);

overall_timer = tic;  % Start timing

for i = 1:stream_count
    [signal, ~, LM_gt, ~, ~, ~, true_labels] = generate_signal_multi_class(Fs, bead_diameters, bead_labels);
    signal = signal - mean(signal);

    % Remove DC offset from signal
    % signal = signal - mean(signal);
    [reco, timing2, LM, RM, Map, Mip]= deriv_method(signal, TH, Fs);
    fprintf('Stream %d: %d detected peaks.\n', i, length(LM));

    if isempty(LM)
        fprintf('  No detected peaks in this stream.\n');
    else
        fprintf('  Detected LM positions (samples):\n');
        disp(LM(1:min(10, length(LM))));
    end
    
    if isempty(LM_gt)
        fprintf('No ground truth labels.\n');
    else
        fprintf('Ground truth LM positions (first 10 samples):\n');
        disp(LM_gt(1:min(10, length(LM_gt))));
    end

    % Filter early peaks
    early_idx = LM > 0 & LM> round(0.1 * Fs); % Exclude peaks in first 10ms
    LM = LM(early_idx);
    RM = RM(early_idx);
    Map = Map(early_idx);
    Mip = Mip(early_idx);

    features = extract_features(LM, RM, Map, Mip, Fs);
    labels = match_detected_to_truth(LM, LM_gt, true_labels);
    fprintf('Stream %d: %d matched labels.\n', i, length(labels));

    % --- Debug plot ---
    debug_plot = true;

    if debug_plot
        figure;
        t = (0:length(signal)-1) / Fs * 1000; % time in ms
        plot(t, signal, 'k'); hold on;
        
        % Plot ground truth
        LM_gt_idx = round(LM_gt);                           % make integers
        LM_gt_idx = LM_gt_idx(LM_gt_idx >= 1 & ...          % keep only legal
                              LM_gt_idx <= numel(signal));
        
        scatter(LM_gt_idx/Fs*1000, signal(LM_gt_idx), ...
                50, 'b', 'filled');
        
        % Plot detected LM (avoid rounding LM)
        LM_ms = round(LM / Fs * 1000);
        scatter(LM_ms, signal(round(LM)), 50, 'r', 'filled');
        
        legend('Signal', 'Ground Truth LM', 'Detected LM');
        xlabel('Time (ms)'); ylabel('Amplitude'); grid on;
        title(sprintf('Stream %d: LM Detected vs. Ground Truth', i)); 
    end

    valid_idx = ~isnan(labels);
    features = features(valid_idx, :);
    labels = labels(valid_idx);
    
    % Confirm alignment
    if size(features, 1) == length(labels)
        all_features = [all_features; features];
        all_labels = [all_labels; labels(:)];  % ensure column vector
    else
        fprintf('Stream %d skipped due to feature-label mismatch.\n', i);
    end

end

if isempty(all_features)
    error('No data collected for training. Check simulation settings.');
end

training_timer = tic;
model = fitctree(all_features, categorical(all_labels));
fprintf('Model training time: %.4f sec\n', toc(training_timer));

save('trained_model.mat', 'model');
fprintf('Classifier retrained and saved as trained_model.mat\n');

fprintf('Total retraining time: %.4f sec\n', toc(overall_timer));
