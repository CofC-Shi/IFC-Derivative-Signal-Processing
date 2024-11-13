function [LL, LM, RM, RR] = find_event_boundaries(peaks, zero_crossings)
    % Initialize event boundary arrays
    LL = zeros(1, length(peaks));
    LM = peaks;
    RM = peaks;
    RR = zeros(1, length(peaks));
    
    for k = 1:length(peaks)
        % Find the closest preceding zero crossing to the left (LL)
        LL_idx = find(zero_crossings < peaks(k), 1, 'last');
        if ~isempty(LL_idx)
            LL(k) = zero_crossings(LL_idx);
        else
            LL(k) = 1; % Set to start if no preceding zero crossing
        end
        
        % Find the closest succeeding zero crossing to the right (RR)
        RR_idx = find(zero_crossings > peaks(k), 1, 'first');
        if ~isempty(RR_idx)
            RR(k) = zero_crossings(RR_idx);
        else
            RR(k) = length(deriv); % Set to end if no succeeding zero crossing
        end
    end
end
