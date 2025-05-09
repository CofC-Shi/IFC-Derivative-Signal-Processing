%% IFC Unsupervised Event Clustering
% Author: Leilei Shi • Rev: May‑2025
clear; close all; clc;

%% ─── USER SETTINGS ────────────────────────────────────────────────
folder_path = "C:\Users\leileishi\OneDrive - College of Charleston\Lab Folder\Leilei\Data\2025\05022025\4um&7umBeadsMix\stream_002";
file_list   = dir(fullfile(folder_path,"*.mat"));     % each file = one stream
TH          = 1e-6;        % derivative threshold for deriv_method
k           = 2;             % #clusters for k‑means  (set to [] for auto‑silhouette)
rng default;                 % reproducibility for k‑means initialisation
%% ─────────────────────────────────────────────────────────────────

% ── collectors
all_amp = [];          % |Map − Mip|
all_tt  = [];          % transit‑time (ms)
all_pkT = [];          % peak times (ms) – for plotting streams

%% ─── PASS 1 : process each .mat file ────────────────────────────
for s = 1:1%numel(file_list)

    disp(['Processing file ', num2str(i), ' of ', num2str(length(file_list)), ': ', file_list(s).name]);
    
    % Load data
    file_path = fullfile(folder_path, file_list(s).name);
    [data, labels, Fs] = LoadData(file_path);

    D = data{2};
    D = D - mean(D);


    % –– detect events
    [~,~,LM,RM,Map,Mip] = deriv_method(D, TH, Fs);
    if isempty(LM),  fprintf('  No events detected.\n');  continue; end

    % –– features
    feats = extract_features(LM, RM, Map, Mip, Fs);  %#ok<NASGU>
    amp   = abs(Map - Mip);               % |Map − Mip| (V)
    tt_ms = (RM - LM)/Fs*1e3;            % transit‑time (ms)

    % –– store
    all_amp = [all_amp; amp(:)];
    all_tt  = [all_tt;  tt_ms(:)];
    all_pkT = [all_pkT; (LM+RM)/2/Fs*1e3];   % mid‑point time (ms)

    % –– quick per‑stream plot (grey trace + event peaks)
    t_ms = (0:numel(D)-1)/Fs*1e3;
    plot_limit = 1000;  % in milliseconds
    figure('Name',['Stream ',num2str(s)]); hold on;
    plot(t_ms, D,'Color',[.7 .7 .7]); xlabel('Time (ms)'); ylabel('V');
    box on                          % enable full box (all sides)
    ax = gca;
    xlim([0 plot_limit]);
    for i = 1:length(LM)
        t1 = LM(i) / Fs * 1000;
        t2 = RM(i) / Fs * 1000;
    
        if t1 > plot_limit, continue; end  % skip if outside plotting window
    
        % Y values for LM and RM
        y1 = D(round(LM(i)));
        y2 = D(round(RM(i)));
    
        % Dashed black line for LM to RM
        plot([t1 t2], [y1 y2], 'Color', [0 0 0 0.3], 'LineWidth', 1.0);
    
        % Get max (positive peak) and min (negative peak) in the range
        range = round(LM(i)):round(RM(i));
        if any(range < 1 | range > numel(D)), continue; end  % safety check
    
        [pos_val, pos_idx] = max(D(range));
        [neg_val, neg_idx] = min(D(range));
        t_pos = (range(1) + pos_idx - 1) / Fs * 1000;
        t_neg = (range(1) + neg_idx - 1) / Fs * 1000;
    
        % Plot peaks
        scatter(t_pos, pos_val, 40, 'r', 'filled');  % red for positive peak
        scatter(t_neg, neg_val, 40, 'b', 'filled');  % blue for negative peak
    end

    title(['Detected events – Stream ',num2str(s)]); grid on; hold off;
end

if isempty(all_amp)
    error('No events found in any file – check paths / thresholds.');
end

%% ─── PASS 2 : k‑means clustering in feature space ───────────────
feat = zscore([all_amp, all_tt]);   % standardise

% pick k automatically by silhouette, if requested
if isempty(k)
    kvals = 2:6;
    sil   = zeros(size(kvals));
    for i = 1:numel(kvals)
        [~,~,~,Dsil] = kmeans(feat,kvals(i),'Replicates',5,'Display','off');
        sil(i) = mean(Dsil);
    end
    [~,best] = max(sil);  k = kvals(best);
    fprintf('Auto‑selected k = %d (silhouette method)\n',k);
end

[idx, Cnorm] = kmeans(feat, k, 'Replicates', 10, 'Display','final');
centroids = Cnorm .* std([all_amp,all_tt]) + mean([all_amp,all_tt]);

%% ─── PLOTS : scatter & histograms ───────────────────────────────
colors = lines(k);
cl_names = "Cluster "+string(1:k);
cluster_cat = categorical(idx, 1:k, cl_names);

% scatter
figure('Name','Transit‑time vs Magnitude (k‑means)');
gscatter(all_amp, all_tt, cluster_cat, colors, 'o', 6);
xlabel('Signal Magnitude (V)');
ylabel('Transit Time  (ms)');
title('Unsupervised Clustering of Events'); grid on;

% histograms
figure('Name','Magnitude histogram (k‑means)'); hold on;
edges = linspace(min(all_amp), max(all_amp), 40);
box on                          % enable full box (all sides)
ax = gca;
for c = 1:k
    mask = idx == c;
    h = histogram(all_amp(mask), edges, 'FaceColor',colors(c,:), ...
                  'FaceAlpha',.45,'DisplayName',cl_names(c));
    [f,xi] = ksdensity(all_amp(mask));
    plot(xi, f*max(h.BinCounts)/max(f), 'Color',colors(c,:), 'LineWidth',1.8);
end
xlabel('Signal Magnitude (V)'); ylabel('# events'); legend; grid on; hold off;

%% ─── PRINT centroids in original units ──────────────────────────
fprintf('\nCluster centroids:\n');
fprintf('  %8s %12s\n','|Map−Mip| (V)','Transit (ms)');
for c = 1:k
    fprintf('C%2d  %10.3g %12.3f\n', c, centroids(c,1), centroids(c,2));
end
