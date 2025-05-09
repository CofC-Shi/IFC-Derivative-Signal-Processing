%% IFC Real-Time Classification System (Modular + Visualization + Retrainable)
% Author: Leilei Shi
% Date: Feb 2025

clear all; close all; clc;

%% Parameters
Fs = 7196*3;           % Hz
duration = 30;                   % seconds per stream
carrier_freq = 1.1e6;            % Hz
TH = 1.5e-6;                      % Derivative threshold
min_P2P_threshold = 0;       % Minimum P2P amplitude
bead_diameters = [6e-6, 5e-6, 4e-6];
bead_labels = [1, 2, 3];        % Class labels

%% Load trained model (must be trained first)
model_path = 'trained_model.mat';
if isfile(model_path)
    load(model_path, 'model');
    disp('Loaded trained model.');
else
    error('Trained model not found. Please train using retrain_classifier().');
end

%% Simulate real-time streaming and classification

% Initialize collection variables
all_transit = [];
all_amplitude = [];
all_true_centers = [];
all_pred_centers = [];
all_predictions = [];
all_true_labels = [];

num_streams = 10;

for stream = 1:num_streams
    fprintf('\n--- Stream %d ---\n', stream);

    %% --- Data Generation ---
    t_start = tic;

    % Simulate new incoming data stream
    [signal, r_clean, LM_gt, RM_gt, Map, Mip, true_labels] = generate_signal_multi_class(Fs, bead_diameters, bead_labels);
    signal = signal - mean(signal);
    t_gen = toc(t_start);
    fprintf('Data generation time: %.4f sec\n', t_gen)

    %% --- Signal Processing (derivative, peak detection, reconstruction) ---
    % Run detection and classification
    t_start = tic;
    [reco, ~, LM, RM, Map, Mip] = deriv_method(signal, TH, Fs);
    t_proc = toc(t_start);
    fprintf('Signal processing time: %.4f sec\n', t_proc);

    %% --- Feature Extraction ---
    t_start = tic;
    features = extract_features(LM, RM, Map, Mip, Fs);
    t_feat = toc(t_start);
    fprintf('Feature extraction time: %.4f sec\n', t_feat);

    %% --- Prediction & Classification ---
    t_start = tic;
    % Prediction
    if ~isempty(features)
        predictions = predict(model, features);
        fprintf('Stream %d: %d particles detected\n', stream, length(predictions));
        disp(['Classes predicted: ', strjoin(string(predictions'), ', ')]);

        transit_time_ms = (RM - LM) / Fs * 1000;
        amplitude = abs(Map - Mip);
        all_transit = [all_transit, transit_time_ms];
        all_amplitude = [all_amplitude, amplitude];

        pred_centers = (LM + RM) / 2 / Fs * 1000;
        true_centers = (LM_gt + RM_gt) / 2 / Fs * 1000;

        all_pred_centers = [all_pred_centers; pred_centers(:)];
        all_true_centers = [all_true_centers; true_centers(:)];
        all_predictions = [all_predictions; predictions];
        all_true_labels = [all_true_labels; categorical(true_labels(:))];
    else
        disp(['Stream ', num2str(stream), ': No events detected.']);
        continue;
    end
    t_pred = toc(t_start);
    fprintf('Classification time: %.4f sec\n', t_pred);

    % --- Live visualization ---
    t_start = tic;
    t_ms = (0:length(signal)-1) / Fs * 1000;
    figure('Name', ['Stream ', num2str(stream)]);
    plot(t_ms, signal, 'Color', [0.6 0.6 0.6]); hold on;
    colors = lines(3);
    
    for c = 1:3
        idx = find(predictions == categorical(c));
        legend_shown = false;  % to prevent repeated legend entries
        
        for j = 1:length(idx)
            k = idx(j);
            if isnan(LM(k)) || isnan(RM(k)) || LM(k) <= 0 || RM(k) > length(signal)
                continue;
            end
            range = round(LM(k)):round(RM(k));
            if isempty(range)
                continue;
            end
            [~, rel_idx] = max(signal(range));
            peak_idx = range(1) + rel_idx - 1;
    
            if ~legend_shown
                scatter((peak_idx - 1)/Fs*1000, signal(peak_idx), ...
                    60, colors(c,:), 'filled', ...
                    'DisplayName', ['Class ', num2str(c)]);
                legend_shown = true;
            else
                scatter((peak_idx - 1)/Fs*1000, signal(peak_idx), ...
                    60, colors(c,:), 'filled', ...
                    'HandleVisibility', 'off');
            end
        end
    end
    
    xlabel('Time (ms)');
    ylabel('Amplitude (V)');
    ylim([-3e-4, 4e-4]);
    title(['Real-Time Classification - Stream ', num2str(stream)]);
    legend('Location','best');
    grid on;
    drawnow;
    t_vis = toc(t_start);
    fprintf('Visualization time: %.4f sec\n', t_vis);

end

%% === Scatter Plot: Transit Time vs. Magnitude (like Figure a)
figure('Name', 'Transit Time vs. Signal Magnitude');
gscatter(all_amplitude, all_transit, all_predictions, 'rgb', 'o', 6);
xlabel('Signal Magnitude |Map + Mip| (V)');
ylabel('Transit Time (ms)');
title('Transit Time vs. Signal Magnitude by Class');
legend('Class 1', 'Class 2', 'Class 3', 'Location', 'northeast');
grid on;

%% === Histogram: Signal Magnitude Distribution (like Figure b)
figure('Name', 'Signal Magnitude Histogram by Class');

edges = linspace(min(all_amplitude), max(all_amplitude), 40);
classes = categories(all_predictions);
colors = lines(length(classes));

hold on;
for i = 1:length(classes)
    class_mask = all_predictions == classes{i};
    h = histogram(all_amplitude(class_mask), edges, ...
        'FaceColor', colors(i,:), ...
        'FaceAlpha', 0.5, ...
        'DisplayName', ['Class ', char(classes{i})]);

    % Optional: KDE curve (smoothed Gaussian fit)
    [f, xi] = ksdensity(all_amplitude(class_mask));
    plot(xi, f * max(h.BinCounts) / max(f), 'Color', colors(i,:), 'LineWidth', 2);
end
xlabel('Signal Magnitude |Map + Mip| (V)');
ylabel('Number of Events');
title('Signal Magnitude Distribution by Class');
legend('Location', 'northeast');
grid on;
hold off;

% === Confusion Matrix ===
%% === Match Predicted Events to Ground Truth Events ===
match_window = 1.0; % ms
matched_preds = [];
matched_truth = [];
used_gt = false(size(all_true_centers));

for i = 1:length(all_pred_centers)
    diffs = abs(all_true_centers - all_pred_centers(i));
    [min_diff, idx] = min(diffs);
    if min_diff <= match_window && ~used_gt(idx)
        matched_preds(end+1) = double(all_predictions(i));
        matched_truth(end+1) = double(all_true_labels(idx));
        used_gt(idx) = true;
    end
end

matched_preds_cat = categorical(matched_preds, bead_labels, "Class " + string(bead_labels));
matched_truth_cat = categorical(matched_truth, bead_labels, "Class " + string(bead_labels));

%% === Confusion Matrix for Matched Events ===
figure;
confusionchart(matched_truth_cat, matched_preds_cat, ...
    'Title', 'Confusion Matrix (Matched Events)', ...
    'RowSummary', 'row-normalized', ...
    'ColumnSummary', 'column-normalized');

%% === Classification Metrics (Precision, Recall, F1) ===
disp('=== Classification Report (Matched Events) ===');
for i = 1:length(bead_labels)
    c = bead_labels(i);
    TP = sum((matched_truth == c) & (matched_preds == c));
    FP = sum((matched_truth ~= c) & (matched_preds == c));
    FN = sum((matched_truth == c) & (matched_preds ~= c));
    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    f1 = 2 * (precision * recall) / (precision + recall + eps);
    fprintf('Class %d: Precision = %.2f, Recall = %.2f, F1-score = %.2f\n', c, precision, recall, f1);
end


