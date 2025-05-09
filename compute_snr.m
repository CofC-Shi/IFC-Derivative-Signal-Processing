%% Event-Based SNR Calculation for IFC Signals
% Author: Leilei Shi (adapted for your data structure)
clear; clc;

% --- Folder Settings ---
noisy_folder = 'mat_data/2.2e-4';
clean_folder = 'mat_data/w1.1e-3_clean';

% --- Load Matching Files ---
noisy_files = dir(fullfile(noisy_folder, '*.mat'));
num_files = length(noisy_files);

% --- Init
SNRs = [];

count = 1;

for i = 1:num_files
    fprintf('Processing %d/%d: %s\n', i, num_files, noisy_files(i).name);

    % Load noisy and clean data
    noisy_path = fullfile(noisy_folder, noisy_files(i).name);
    clean_path = fullfile(clean_folder, noisy_files(i).name);

    noisy_data = LoadData(noisy_path);
    clean_data = LoadData(clean_path);

    % Extract channel 2
    sig_noisy = noisy_data{2};
    sig_noisy = sig_noisy - mean(sig_noisy);
    sig_clean = clean_data{2};
    sig_clean = sig_clean - mean(sig_clean);

    %pull event windows
    events = detect_ev(sig_clean,0);
    av_trans = mean(events(2:2:end) - events(1:2:end-1));
    %Calculate SNR inside event windows
    SNRs(i) = CalcSNR(sig_noisy,sig_clean,events);
    disp([num2str(count),'/',num2str(num_files)])
    count = count + 1;

end

% --- Report ---
SNR_total = mean(SNRs);
fprintf('Average Event-Based SNR: %.2f dB\n', SNR_total);


function SNRs = CalcSNR(sig,clean,events)

    win = length(events)/2;
    noise = sig - clean;
    p_sp = [];
    p_np = [];
    
    i2 = 0;
    for k = 1:win
        i1 = 1 + i2;
        i2 = 2*k;
        index = events(i1):events(i2);
        p_sp(k) = rms(clean(index))^2;
        p_np(k) = rms(noise(index))^2;
    end
    
    p_s = mean(p_sp);
    p_n = mean(p_np);
    
    SNRs = 10*log10(p_s/p_n);

end