% Title: Peak Detection Notch
% Author: Leilei Shi
% Date: 11-08-24
%
% Description: Loads, processes, and detects events in raw datastreams 
%              using a simple notch filtering, peak detection, and thresholding method.
%              After processing, extract positive peak, negative peak, and the
%              time peak to peak.

clc;
close all;
clear all;

% Folder path
folder_path = "Input Path";

% List all .mat files in the folder
file_list = dir(fullfile(folder_path, "*.mat"));

% Set sensitivity for peak detection
sensitivity = 0.2;
noise_threshold_multiplier = 4;  % Minimum peak height as a multiple of noise level

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
        % % Display progress for each channel
        % disp(['  Processing channel ', num2str(channel), ' of ', num2str(length(data))]);
        d = data{channel};
        
        % Filter and detrend data
        processed_data = DataProcess(d, Fs);

        % Estimate noise level as the standard deviation of the signal
        noise_level = std(processed_data);
        min_peak_height = noise_threshold_multiplier * noise_level;
        
        % Find Positive Peaks
        data_max = max(processed_data);
        [pks_pos, locs_pos] = findpeaks(processed_data, Fs, "MinPeakHeight", min_peak_height);

        % Filter positive peaks to include only those after 1 second
        locs_pos = locs_pos(locs_pos > 1); % Only include peaks after 1 second
        pks_pos = pks_pos(locs_pos > 1); % Corresponding positive peak values
        
        % Find Negative Peaks (invert data)
        [pks_neg, locs_neg] = findpeaks(-processed_data, Fs, "MinPeakHeight", min_peak_height);
        pks_neg = -pks_neg;  % Revert to original negative values

        % % Filter negative peaks to include only those after 1 second
        % locs_neg = locs_neg(locs_neg > 1); % Only include peaks after 1 second
        % pks_neg = pks_neg(locs_neg > 1); % Corresponding positive peak values

        % Skip processing if detected peaks do not stand out from noise
        if isempty(pks_pos) || isempty(pks_neg)
            disp(['No significant signal detected in ', file_list(i).name, ', Channel: ', num2str(channel)]);
            continue;
        end
        
        % Pair positive and negative peaks within 1-second window
        % Pair positive and negative peaks within a 0.5-second window
        valid_locs_pos = [];
        valid_pks_pos = [];
        valid_locs_neg = [];
        valid_pks_neg = [];
        peak_to_peak_times = [];
        
        % Loop through each negative peak
        for j = 1:length(locs_neg)
            % Find the closest positive peak that occurs after the current negative peak and within 0.5 seconds
            subsequent_pos_indices = find(locs_pos > locs_neg(j) & (locs_pos - locs_neg(j)) <= 0.5);
            
            % Check if there's a valid positive peak within the window
            if isempty(subsequent_pos_indices)
                % If there is no valid positive peak within 0.5 seconds, skip this negative peak
                continue; 
            end
            
            % Select the nearest positive peak within the window
            nearest_pos_index = subsequent_pos_indices(1);
            
            % Calculate the peak-to-peak time
            peak_to_peak_time = abs(locs_neg(j) - locs_pos(nearest_pos_index)) * 1000; % Convert to ms
            
            % Store the valid negative and positive peak pair and peak-to-peak time
            valid_locs_neg = [valid_locs_neg; locs_neg(j)];
            valid_pks_neg = [valid_pks_neg; pks_neg(j)];
            valid_locs_pos = [valid_locs_pos; locs_pos(nearest_pos_index)];
            valid_pks_pos = [valid_pks_pos; pks_pos(nearest_pos_index)];
            peak_to_peak_times = [peak_to_peak_times; peak_to_peak_time];
        end

        % Check if any valid pairs were found
        if isempty(valid_locs_pos) || isempty(valid_locs_neg)
            disp(['No valid peak pairs found within 1 second in ', file_list(i).name, ', Channel: ', num2str(channel)]);
            continue;  % Skip further processing for this file/channel if no valid pairs
        end

        % Ensure both positive and negative peaks have the same length by trimming
        min_length = min([length(valid_locs_pos), length(valid_locs_neg)]);
        valid_locs_pos = valid_locs_pos(1:min_length);
        valid_pks_pos = valid_pks_pos(1:min_length);
        valid_locs_neg = valid_locs_neg(1:min_length);
        valid_pks_neg = valid_pks_neg(1:min_length);
        peak_to_peak_times = peak_to_peak_times(1:min_length);

        % **Condition to skip saving if too many peaks are detected (noise threshold)**
        if length(valid_pks_pos) + length(valid_pks_neg) > 30
            disp(['Excessive noise detected in ', file_list(i).name, ', Channel: ', num2str(channel), ' - Skipping file.']);
            continue;  % Skip saving if total peaks exceed threshold
        end
         % End timing for this file
        elapsed_time = toc;
        
        % Add the elapsed time to the total time
        total_time = total_time + elapsed_time;
        
        % Display the processing time for this file
        disp(['Processing time for ', file_list(i).name, ': ', num2str(elapsed_time), ' seconds']);
        
        % Create table with data for this channel
        channel_data = table(peak_to_peak_times', valid_pks_pos', valid_pks_neg', ...
                             'VariableNames', {'Time_ms', 'PosPeak', 'NegPeak'});

        % Define unique output file name for each channel
        output_file = fullfile(folder_path, ...
            sprintf('%s_Channel%d_results_4um.csv', file_list(i).name(1:end-4), channel));

        % Save channel data to CSV file
        writetable(channel_data, output_file);
        disp(['Results saved to ', output_file]);

        % Plot the original and processed data with detected peaks
        figure;
        time = (0:length(d)-1) / Fs * 1000;  % Time in ms

        % Top subplot: Original signal
        subplot(2, 1, 1);
        plot(time, d, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  % Original signal in gray
        legend('Original Data')
        xlim([0, 1e4]);
        hold on;

        % Plot the processed data with detected peaks
        subplot(2, 1, 2);
        time = (0:length(processed_data)-1) / Fs * 1000;  % Time in ms
        plot(time, processed_data);
        hold on;
        scatter(valid_locs_pos*1000, valid_pks_pos, 'r', 'filled');  % Positive peaks in red
        scatter(valid_locs_neg*1000, valid_pks_neg, 'b', 'filled');  % Negative peaks in blue
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