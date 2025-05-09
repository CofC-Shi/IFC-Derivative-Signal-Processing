function labels = match_detected_to_truth(LM_detected, LM_true, true_labels)
% match_detected_to_truth - Matches detected events to true event locations based on nearest neighbor matching
%
% Inputs:
%   LM_detected - Detected left midpoints (from peak detection)
%   LM_true - Ground truth left midpoints (from simulation)
%   true_labels - Class labels of ground truth events
%
% Output:
%   labels - Matched labels for detected events

    tolerance = 40;  % samples allowed (adjustable)
    labels = NaN(size(LM_detected));  % initialize with NaNs
    assigned = false(size(LM_true));  % track which ground truth events have been matched

    for i = 1:length(LM_detected)
        diffs = abs(LM_true - LM_detected(i));
        diffs(assigned) = inf;  % ignore already assigned ground truth
        [min_dist, idx] = min(diffs);
        if min_dist <= tolerance
            labels(i) = true_labels(idx);
            assigned(idx) = true;
        end
    end

    % Remove unmatched
    valid_idx = ~isnan(labels);
    labels = labels(valid_idx);
end
