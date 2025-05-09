function [data_out, freq_labels, sample_freq] = LoadData(data_address)

% DESCRIPTION
%
%   Loads and formats data from .mat file.
%
%   Author - Charlie Jindrich
%   Date - 05/07/24

% INPUTS
%
%   data_address - string containing the file path of data

% OUTPUTS
%  
%   data_out - matrix with each column representing impedance data of 
%              a frequency channel.
%
%   freq_labels - a cell array of frequencies for each channel.
%
%   sample_freq - Sampling frequency of data.

% Load data structure
data_struct = load(data_address);

% Extract Sample Freq
sample_freq = data_struct.dev18244.demods(1).rate.value;

% Find number of frequency channels
num_elements = numel(data_struct.dev18244.demods);

% Preallocate array for data
data_out = cell(1, num_elements);
freq_labels = cell(1, num_elements);

len = length(data_struct.dev18244.demods(1).sample.x);

% Extract each channel
for i = 1:num_elements

    x = data_struct.dev18244.demods(i).sample.x;
    y = data_struct.dev18244.demods(i).sample.y;
    r = sqrt((x.^2) + (y.^2));

    % Make sure all channels are the same length
    if length(r) ~= len
        if length(r) < len
            while length(r) < len
                r = horzcat(r, mean(r));
            end
        else
            r = r(1:len);
        end
    end

    data_out{i} = r;
    freq_labels{i} = data_struct.dev18244.demods(i).freq.value;

end
end