%% IFC Synthetic Signal Generator - Single-Class Version
% Author: Simplified for Leilei Shi (One class per file)
clear; clc;

% Parameters
num_streams = 30;              % number of .mat files to generate
Fs = 7196*5;                    % sampling rate
output_folder_clean = 'mat_data/fast_20_clean/';
output_folder_noisy = 'mat_data/fast_20/';
mkdir(output_folder_clean);
mkdir(output_folder_noisy);

% Bead class: adjust for each run
bead_diameter = 7.32e-6;
bead_label = 1;

%% Generation loop
disp("Generating .mat files for LoadData.m...");

for k = 1:num_streams
    % Generate signal with noise and clean reference
    [r_noisy, r_clean, LM_gt, RM_gt, Map, Mip, ~] = generate_signal_multi_class(Fs, bead_diameter, bead_label);

    % Simulate approximate x/y assuming zero phase
    x_noisy = r_noisy;
    y_noisy = zeros(size(r_noisy));
    x_clean = r_clean;
    y_clean = zeros(size(r_clean));

    % Create noisy structure
    dev18244 = struct();
    for ch = 1:2
        dev18244.demods(ch).sample.x = x_noisy;
        dev18244.demods(ch).sample.y = y_noisy;
        dev18244.demods(ch).rate.value = Fs;
        dev18244.demods(ch).freq.value = 0.5e6 * ch;
    end
    save(fullfile(output_folder_noisy, sprintf('stream_simulated_%02d.mat', k))), 'dev18244';

    % Create clean structure
    dev18244= struct(); % reuse structure layout
    for ch = 1:2
        dev18244.demods(ch).sample.x = x_clean;
        dev18244.demods(ch).sample.y = y_clean;
        dev18244.demods(ch).rate.value = Fs;
        dev18244.demods(ch).freq.value = 0.5e6 * ch;
    end
    % save(fullfile(output_folder_clean, sprintf('stream_simulated_%02d.mat', k))), 'dev18244';
    save(fullfile(output_folder_clean, sprintf('stream_simulated_%02d.mat', k)), ...
     'dev18244', 'LM_gt', 'RM_gt', 'Map', 'Mip');


    disp(['Stream ', num2str(k), ' saved.']);
end

disp("All done!");
