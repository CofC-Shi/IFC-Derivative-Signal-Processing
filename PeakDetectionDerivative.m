% Title: Clean Peak Detection with Derivative Method
% Author: Leilei Shi
% Date: 11-07-2024
%
% Description:
% Processes signals using derivative method.
% Extracts positive/negative peaks from LM, RM directly.
% Builds clean table and generates organized plots.

clear all;
close all;
clc;

% --- Settings --1
folder_path = "mat_data/fast_20";
clean_folder = folder_path + "_clean";
clean = true;
file_list = dir(fullfile(folder_path, "*.mat"));

% folder_path = "C:\Users\leileishi\OneDrive - College of Charleston\Lab Folder\Leilei\Data\2025\05022025\4um&7umBeadsMix\stream_002";
% clean = false;
% file_list = dir(fullfile(folder_path, "*.mat"));

TH = 1.5e-6; % Threshold for derivative-based peak detection
min_P2P_threshold = 0%1e-5; % The minimal peak-to-peak value

% Initialize total time
total_time = 0;

% --- Process Each File ---
for i = 1:1 %length(file_list) % Change 1 to full batch when needed
    tic; % Start timer
    
    disp(['Processing file ', num2str(i), ' of ', num2str(length(file_list)), ': ', file_list(i).name]);
    
    % Load data
    file_path = fullfile(folder_path, file_list(i).name);
    [data, labels, Fs] = LoadData(file_path);

    d = data{2};
    d = d - mean(d);

    % --- Load Clean Signal ---
    if clean
        % Extract the name part (e.g., '1e-5') and build the clean path
        clean_path = fullfile(clean_folder, file_list(i).name);
    
        disp(['Clean file path: ', clean_path]);  % <-- Print the clean path
        clean_data = LoadData(clean_path);
        sig_clean = clean_data{2};
        sig_clean = sig_clean - mean(sig_clean);
    end

    % --- Apply derivative-based method ---
    [reco, timing2, LM, RM, Map, Mip] = deriv_method(d, TH, Fs);

    % --- Extract Positive/Negative Peaks ---

    % Convert LM and RM to milliseconds
    Time_Positive_ms = (LM(:) - 1) / Fs * 1000; % column vector
    Time_Negative_ms = (RM(:) - 1) / Fs * 1000;
    Amplitude_Positive = Map(:);
    Amplitude_Negative = Mip(:);
    PeakToPeakAmplitude = Amplitude_Positive - Amplitude_Negative;

    % Calculate Peak-to-Peak Times
    PeakToPeakTime_ms = (RM(:) - LM(:)) / Fs * 1000;

    % valid_idx = (Amplitude_Positive >= 0.5 * min_P2P_threshold) & ...
    %     (Amplitude_Negative <= -0.5 * min_P2P_threshold);
    % 
    % % Keep only valid events
    % Time_Positive_ms = Time_Positive_ms(valid_idx);
    % Amplitude_Positive = Amplitude_Positive(valid_idx);
    % Time_Negative_ms = Time_Negative_ms(valid_idx);
    % Amplitude_Negative = Amplitude_Negative(valid_idx);
    % PeakToPeakAmplitude = PeakToPeakAmplitude(valid_idx);
    % PeakToPeakTime_ms = PeakToPeakTime_ms(valid_idx);
    num_events = length(Time_Positive_ms); % Update number of events

    % % --- Check if any valid peaks ---
    % if isempty(Time_Positive_ms)
    %     disp(['No valid peaks found in ', file_list(i).name, ', skipping...']);
    %     continue;
    % end


    % --- Build Table ---
    channel_data = table(Time_Positive_ms, ...
                         Amplitude_Positive, ...
                         Time_Negative_ms, ...
                         Amplitude_Negative, ...
                         PeakToPeakAmplitude, ...
                         PeakToPeakTime_ms, ...
                         'VariableNames', {'Time_Positive_ms', 'Amplitude_Positive', ...
                                           'Time_Negative_ms', 'Amplitude_Negative', ...
                                           'PeakToPeakAmplitude','PeakToPeakTime_ms'});

    % --- Save Table ---
    output_file = fullfile(folder_path, ...
        sprintf('%s_Channel%d_results_7um_derivative.csv', file_list(i).name(1:end-4), 2));
    writetable(channel_data, output_file);
    disp(['Results saved to ', output_file]);

  % --- Plotting ---
    
  plotting = true;
  if plotting

     %% 3 rows plot: original, derivative1, derivative2, reconstructed
    figure;
    tiledlayout(3, 1);  
    t.TileSpacing = 'tight';     % No vertical spacing
    t.Padding = 'compact';       % Minimal outer margins
    
    % --- Top subplot: Original signal ---
    nexttile;
    time_original = (0:length(d)-1) / Fs * 1000;  % ms
    
    plot(time_original, d, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    if clean
        hold on;
        plot(time_original, sig_clean(1:length(d)), 'g');  % Green = clean signal
    end
    % xlim([120 1000]);
    % ylim([-4e-4 4e-4 ]);
    xlabel('Time (ms)');
    ylabel('Amplitude (V)');
    title('Original Signal');
    grid on;
    
    % --- 2nd subplot: First derivative ---
    nexttile;
    deriv = BDC(d);
    time_vector = (0:length(deriv)-1) / Fs * 1000;  % ms

    plot(time_vector, deriv, 'b');
    % xlim([120 1000]);
    xlabel('Time (ms)');
    ylabel('dV/dt (V/s)');
    title('Derivative of Signal');
    grid on;
    
    % --- Bottom subplot: Reconstructed signal + events ---
    nexttile;
    time_reco = (0:length(reco)-1) / Fs * 1000;  % ms
    p1 = plot(time_reco, reco, 'LineWidth', 1.2);
    set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    % xlim([120 1000]);
    % ylim([-4e-4 4e-4 ]);
    hold on;
    
    scatter(Time_Positive_ms, Amplitude_Positive, 'r', 'filled');
    scatter(Time_Negative_ms, Amplitude_Negative, 'b', 'filled');
    
    for k = 1:num_events
        plot([Time_Positive_ms(k), Time_Negative_ms(k)], ...
             [Amplitude_Positive(k), Amplitude_Negative(k)], ...
             '--k', 'LineWidth', 1.5);
    end
    
    xlabel('Time (ms)');
    ylabel('Amplitude (V)');
    legend('Positive Peak', 'Negative Peak', 'Peak-to-Peak Connection');
    title('Reconstructed Signal');
    grid on;
    hold off;

   %% --- Frequency Spectrum Plot ---
    n = length(d);
    f = Fs*(0:(n/2))/n;
    
    % FFT for Original Signal
    Y_orig = fft(d);
    P2_orig = abs(Y_orig/n);
    P1_orig = P2_orig(1:n/2+1);
    P1_orig(2:end-1) = 2*P1_orig(2:end-1);
    
    % FFT for Reconstructed Signal
    Y_reco = fft(reco);
    P2_reco = abs(Y_reco/n);
    P1_reco = P2_reco(1:n/2+1);
    P1_reco(2:end-1) = 2*P1_reco(2:end-1);
    
    % Avoid log(0)
    P1_orig(P1_orig <= 0) = 1e-12;
    P1_reco(P1_reco <= 0) = 1e-12;
    f(f == 0) = 1e-12;
    
    % --- Combined Frequency Spectrum Plot with Separator ---

    figure;
    t = tiledlayout(2, 1, 'TileSpacing', 'none', 'Padding', 'none');
    
    % Top: Original
    nexttile;
    plot(f, P1_orig, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    ylabel('Amplitude');
    xlim([10 500]);
    % ylim([0, max(P1_orig)*1.1]);
    ylim([0, 0.3e-4]);
    grid on;
    set(gca, 'XTickLabel', [])
    
    % Bottom: Reconstructed
    nexttile;
    plot(f, P1_reco*1e5, 'b', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    xlim([10 500]);
    ylim([0, 3]);
    grid on;

    
    %% --- Spectrograms ---
    figure;
    tiledlayout(2,1);
    
    nexttile;
    spectrogram(d, 512, 500, 1024, Fs, 'yaxis');
    title('Spectrogram - Original Signal');
    ylim([0 3.5]);
    colorbar;
    
    nexttile;
    spectrogram(reco, 512, 500, 1024, Fs, 'yaxis');
    title('Spectrogram - Reconstructed Signal');
    ylim([0 3.5]);
    colorbar;
    
    %% --- Noise Energy Distribution (Smooth Line Plot) ---

    % Define 100 evenly spaced frequency bands from 0 to 500 Hz
    num_bands = 50;
    edges = linspace(0, 1000, num_bands + 1);   % 101 edges define 100 bands
    bands = [edges(1:end-1)', edges(2:end)'];  % 100 rows of [start, end]
    f_centers = mean(bands, 2);                % Center frequencies
    
    % Compute energy in each band
    get_band_energy = @(f, P, band) sum(P(f >= band(1) & f < band(2)).^2);
    energy_orig = arrayfun(@(i) get_band_energy(f, P1_orig, bands(i,:)), 1:num_bands);
    energy_reco = arrayfun(@(i) get_band_energy(f, P1_reco, bands(i,:)), 1:num_bands);
    
    % Normalize energy values
    energy_orig = energy_orig / sum(energy_orig);
    energy_reco = energy_reco / sum(energy_reco);
    
    % Smooth interpolation for curves
    f_dense = linspace(min(f_centers), max(f_centers), 1000);
    energy_orig_interp = interp1(f_centers, energy_orig, f_dense, 'pchip');
    energy_reco_interp = interp1(f_centers, energy_reco, f_dense, 'pchip');
    
    % Plot
    figure;
    hold on;

    bar(f_centers, energy_orig, 0.8, 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'FaceColor', 'b');
    bar(f_centers, energy_reco, 0.8, 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'FaceColor', 'r');


    % Overlay smooth lines
    plot(f_dense, energy_orig_interp, 'b-', 'LineWidth', 1);
    plot(f_dense, energy_reco_interp, 'r-', 'LineWidth', 1);

    xlabel('Frequency (Hz)');
    ylabel('Normalized Energy');
    title('Noise Energy Distribution (0–500 Hz)');
    legend('Original (hist)', 'Reconstructed (hist)', 'Original (curve)', 'Reconstructed (curve)');
    grid on;
    xlim([0 500]);

    % % % --- Dual-Y Plot: Histogram + Smooth Curves ---
    % 
    % figure;
    % 
    % % === Histogram Bars ===
    % yyaxis left
    % bar(f_centers, energy_orig, 0.8, ...
    %     'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [0.4 0.4 1]);  % light blue
    % ylabel('Original Normalized Energy');
    % ylim([0 max(energy_orig)*1.2]);
    % 
    % yyaxis right
    % bar(f_centers, energy_reco, 0.8, ...
    %     'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [1 0.4 0.4]);  % light red
    % ylabel('Reconstructed Normalized Energy');
    % ylim([0 max(energy_reco)*1.2]);
    % 
    % % === Overlay smooth lines ===
    % hold on;
    % yyaxis left;
    % plot(f_dense, energy_orig_interp, 'b-', 'LineWidth', 1);
    % 
    % yyaxis right;
    % plot(f_dense, energy_reco_interp, 'r-', 'LineWidth', 1);
    % 
    % % === Axes Labels and Aesthetics ===
    % xlabel('Frequency (Hz)');
    % title('Noise Energy Distribution (0–500 Hz)');
    % legend({'Original (hist)', 'Reconstructed (hist)', ...
    %         'Original (curve)', 'Reconstructed (curve)'}, ...
    %         'Location', 'northoutside', 'Orientation', 'horizontal');
    % xlim([0 500]);
    % grid on;
  end

    % Timer
    elapsed_time = toc;
    total_time = total_time + elapsed_time;
    disp(['Processing time for ', file_list(i).name, ': ', num2str(elapsed_time), ' seconds']);
end

% Average time
average_time = total_time / length(file_list);
disp(['Average processing time per file: ', num2str(average_time), ' seconds']);
