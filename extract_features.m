function features = extract_features(LM, RM, Map, Mip, Fs)
% extract_features - Extract peak-to-peak amplitude and transit time from detected peaks
%
% Inputs:
%   LM - Left midpoints of peaks
%   RM - Right midpoints of peaks
%   Map - Positive peak amplitudes
%   Mip - Negative peak amplitudes
%   Fs - Sampling frequency
%
% Outputs:
%   features - [PeakToPeakAmplitude, PeakToPeakTime_ms]

    if isempty(LM)
        features = [];
        return;
    end

    PeakToPeakAmplitude = Map(:) - Mip(:);  % Difference between positive and negative peaks
    PeakToPeakTime_ms = (RM(:) - LM(:)) / Fs * 1000;  % Peak-to-peak time in milliseconds

    features = [PeakToPeakAmplitude, PeakToPeakTime_ms];
end
