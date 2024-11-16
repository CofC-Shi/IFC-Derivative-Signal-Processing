% Title: Peak Detection Notch
% Author: Leilei Shi
% Date: 11-15-24
%
% Description: Loads, processes, and detects events in raw datastreams 
%              using a simple notch filtering, peak detection, and thresholding method.
%              After processing, extract positive and negative peaks separately
%              along with their widths and times.


clc;
close all;
clear all;

% Folder path
folder_path = "Input Path";

% List all .mat files in the folder
file_list = dir(fullfile(folder_path, "*.mat"));

% Set threshold for peak detection
noise_threshold = 0.5e-4;  % Minimum peak height as a multiple of noise level
% Define sampling parameters
sampling_frequency = 7196; % in Hz

% Initialize a variable to store total processing time
total_time = 0;

% Process each file
for i = 1:length(file_list)
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
        
        % Find Positive Peaks
        [pks_pos, locs_pos, w_pos, p_pos] = findpeaks(processed_data, Fs, "MinPeakHeight", noise_threshold, 'MinPeakDistance',0.1);
        
        % Find Negative Peaks (invert data)
        [pks_neg, locs_neg, w_neg, p_neg] = findpeaks(-processed_data, Fs, "MinPeakHeight", noise_threshold, 'MinPeakDistance',0.1);
        pks_neg = -pks_neg;  % Revert to original negative values

       % Filter peaks with prominence
        valid_pos_idx = p_pos > 1.5e-4;
        valid_neg_idx = p_neg > 1.5e-4;
        
        % Filter positive peaks
        pks_pos_filtered = pks_pos(valid_pos_idx);
        locs_pos_filtered = locs_pos(valid_pos_idx);
        w_pos_filtered = w_pos(valid_pos_idx);
        p_pos_filtered = p_pos(valid_pos_idx);
        
        % Filter negative peaks
        pks_neg_filtered = pks_neg(valid_neg_idx);
        locs_neg_filtered = locs_neg(valid_neg_idx);
        w_neg_filtered = w_neg(valid_neg_idx);
        p_neg_filtered = p_neg(valid_neg_idx);

        % Convert times and widths to milliseconds
        locs_pos_ms = locs_pos_filtered / sampling_frequency * 1000;
        locs_neg_ms = locs_neg_filtered / sampling_frequency * 1000;
        widths_pos_ms = w_pos_filtered / sampling_frequency * 1000;
        widths_neg_ms = w_neg_filtered / sampling_frequency * 1000;

        % stack positive and negative signals together
        locs_ms = [locs_pos_ms, locs_neg_ms];
        pks_filtered = [pks_pos_filtered, -pks_neg_filtered];
        widths_ms = [widths_pos_ms, widths_neg_ms];
        p_filtered = [p_pos_filtered, p_neg_filtered];

        
        % End timing for this file
        elapsed_time = toc;
        
        % Add the elapsed time to the total time
        total_time = total_time + elapsed_time;
        
        % Display the processing time for this file
        disp(['Processing time for ', file_list(i).name, ': ', num2str(elapsed_time), ' seconds']);

       % Create the final table with separate columns for each parameter
        channel_data = table(locs_ms, pks_filtered, widths_ms, p_filtered, ...
                             'VariableNames', {'Time_ms', 'Amplitude', 'Width', 'Prominence'});

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
        % xlim([0, 1e4]);
        hold on;

        % Plot the processed data with detected peaks
        subplot(2, 1, 2);
        time = (0:length(processed_data)-1) / Fs * 1000;  % Time in ms
        plot(time, processed_data);
        hold on;
        scatter(locs_pos_filtered*1000, pks_pos_filtered, 'r', 'filled');  % Positive peaks in red
        scatter(locs_neg_filtered*1000, pks_neg_filtered, 'b', 'filled');  % Negative peaks in blue
        % title(['Detected Peaks - File: ', file_list(i).name, ', Channel: ', num2str(channel)]);
        xlabel('Time (ms)');
        ylabel('Amplitude (V)');
        legend('Processed Data', 'Positive Peaks', 'Negative Peaks');
        % xlim([0, 1e4]);
        % ylim([-1.5e-4, 3e-4]);
        hold off;
    end
    % Calculate and display the average processing time per file
    average_time = total_time / length(file_list);
    disp(['Average processing time per file: ', num2str(average_time), ' seconds']);
end
