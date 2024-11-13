% Title: Peak Detection with Derivative Method
% Author: Leilei Shi
% Date: 11-07-2024
%
% Description: Loads, processes, and detects positive and negative peaks in raw datastreams
%              using a derivative-based method, and calculates peak-to-peak times.

clear all;
close all;
clc;

% Folder path
folder_path = "Input Folder";

% List all .mat files in the folder
file_list = dir(fullfile(folder_path, "*.mat"));

% Initialize a variable to store total processing time
total_time = 0;

% Process each file
for i = 10:20%length(file_list)
    % Start timing for this file
    tic;
    % Display progress for each file
    disp(['Processing file ', num2str(i), ' of ', num2str(length(file_list)), ': ', file_list(i).name]);
    
    % Load data file
    file_path = fullfile(folder_path, file_list(i).name);
    [data, labels, Fs] = LoadData(file_path);
    
    % Process each channel in the data
    for channel = 2:2
        d_detrended = -data{channel};

        % Estimate noise level as the standard deviation of the processed signal
        % d_detrended = detrend(d);
        % noise_level = std(d_detrended);
        % % disp(['Noise level: ', num2str(noise_level)]);

        % Set threshold based on noise level
        TH = 0.3

        % Define sampling parameters
        sampling_frequency = 7196; % in Hz
        
        % Apply the derivative-based method to the data
        [processed_data, der_t, LM, RM, Map, Mip] = deriv_method(d_detrended, TH);
        locs_neg = RM;
        pks_neg = Mip;
        locs_pos = LM;
        pks_pos = Map;

        % Convert indices to time
        time_pos = locs_pos / sampling_frequency; % Positive peak times in seconds
        time_neg = locs_neg / sampling_frequency; % Negative peak times in seconds

        % Filter out pairs where the negative peak is not less than 0
        valid_pairs = pks_neg < 0;  % Logical array where true indicates valid positive peaks
        time_pos = time_pos(valid_pairs);
        time_neg = time_neg(valid_pairs);
        pks_pos = pks_pos(valid_pairs);
        pks_neg = pks_neg(valid_pairs);

        % **Condition to skip saving if too many peaks are detected (noise threshold)**
        if length(pks_pos) + length(pks_neg) > 50
            disp(['Excessive noise detected in ', file_list(i).name, ', Channel: ', num2str(channel), ' - Skipping file.']);
            continue;  % Skip saving if total peaks exceed threshold
        end
        
        % Ensure both locs_pos and locs_neg are matched in pairs
        num_peaks = min(length(time_pos), length(time_neg));
        time_pos = time_pos(1:num_peaks);
        time_neg = time_neg(1:num_peaks);
        pks_pos = pks_pos(1:num_peaks);
        pks_neg = pks_neg(1:num_peaks);

        % Calculate peak-to-peak time
        peak_to_peak_time = -1*(time_pos(1:num_peaks) - time_neg(1:num_peaks)) * 1000; % Convert to ms

        % End timing for this file
        elapsed_time = toc;
        
        % Add the elapsed time to the total time
        total_time = total_time + elapsed_time;
        
        % Display the processing time for this file
        disp(['Processing time for ', file_list(i).name, ': ', num2str(elapsed_time), ' seconds']);

        % Create table with data for this channel
        channel_data = table(peak_to_peak_time, pks_pos, pks_neg, ...
                             'VariableNames', {'Time_ms', 'PosPeak', 'NegPeak'});

        % Define unique output file name for each channel
        output_file = fullfile(folder_path, ...
            sprintf('%s_Channel%d_results_4um_derivative.csv', file_list(i).name(1:end-4), channel));

        % Save channel data to CSV file
        writetable(channel_data, output_file);
        disp(['Results saved to ', output_file]);

        Plot the original and processed data with detected peaks
        figure;
        time = (0:length(d_detrended)-1) / Fs * 1000;  % Time in ms

        % Top subplot: Original signal
        subplot(2, 1, 1);
        plot(time, -d_detrended, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  % Original signal in gray
        legend('Original Data')
        xlim([0, 1e4]);
        hold on;

        % Plot the processed data with detected peaks
        subplot(2, 1, 2);
        time = (0:length(processed_data)-1) / Fs * 1000;  % Time in ms
        plot(time, processed_data);
        hold on;
        scatter(time_pos * 1000, pks_pos, 'r', 'filled');  % Positive peaks in red
        scatter(time_neg * 1000, pks_neg, 'b', 'filled');  % Negative peaks in blue
        % title(['Detected Peaks - File: ', file_list(i).name, ', Channel: ', num2str(channel)]);
        xlabel('Time (ms)');
        ylabel('Amplitude (V)');
        legend('Processed Data', 'Positive Peaks', 'Negative Peaks');
        xlim([0, 1e4]);
        ylim([-1.5e-4, 3e-4]);
        hold off;
    end
    % Calculate and display the average processing time per file
    average_time = total_time / length(file_list);
    disp(['Average processing time per file: ', num2str(average_time), ' seconds']);
end
