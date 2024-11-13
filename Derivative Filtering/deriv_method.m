% Derivative Method for Signal Processing
% This function processes the input signal using a derivative-based method to identify events,
% reconstructs the signal, and measures the execution time.

function [reco, timing2, LM, RM, Map, Mip] = deriv_method(sig, TH)

% Start timing the execution
t_start = clock;

% % Gaussian kernel
% sigma = 1.5; % Standard deviation for Gaussian kernel
% kernel_size = 7; % Size of the Gaussian kernel (choose an odd number)
% 
% % Create Gaussian kernel
% x = -floor(kernel_size/2):floor(kernel_size/2);
% gaussian_kernel = exp(-x.^2 / (2 * sigma^2));
% gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize
% 
% % Take the derivative of the Gaussian kernel
% d_gaussian_kernel = gradient(gaussian_kernel);
% 
% % Convolve the signal with the derivative of the Gaussian kernel
% deriv = convn(sig, d_gaussian_kernel, 'same');

% % Calculate the derivative of the signal using a custom BDC function
deriv = BDC(sig);

% Define sampling parameters
num_samples = length(sig);
sampling_frequency = 7196; % in Hz
time_vector = (0:num_samples-1) / sampling_frequency; % Time vector in seconds

% % Plot the Derivative of the Signal
% plot(time_vector, deriv * 1000); % Convert to mV/s for visualization
% title('Derivative of Signal');
% xlabel('Time (s)');
% ylabel('Derivative (mV/s)');

% Detect Peaks in the Derivative Signal
[px, py] = NP(TH, deriv);

% Find Zero Crossings for Labeling Events
[LL, LM, RM, RR] = LaR(px, deriv);

% Calculate Peak Values from Detected Events
[Map, Mip] = calcpeak(LL, LM, RM, deriv);

% Reconstruct the Signal using a Bi-Gaussian Model
reco = bigaus_recon(LM, RM, Map, Mip, deriv);


% Measure the Execution Time in Milliseconds
timing2 = round(1000 * etime(clock, t_start));

end


