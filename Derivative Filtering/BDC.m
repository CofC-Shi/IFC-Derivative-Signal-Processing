function diff = BDC(signal)
    sig_len = length(signal);

    % Preallocate the diff array
    diff = zeros(1, sig_len);

    % Forward difference for the first point
    diff(1) = signal(2) - signal(1);

    % Central differences for the interior points
    diff(2:sig_len-1) = (signal(3:sig_len) - signal(1:sig_len-2)) / 2;

    % Backward difference for the last point
    diff(sig_len) = signal(sig_len) - signal(sig_len - 1);
end
