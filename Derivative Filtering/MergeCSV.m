% Folder path where CSV files are located
folder_path = "Output Path";

% List all CSV files in the folder
file_list = dir(fullfile(folder_path, "*_results_4um_derivative.csv"));

% Initialize an empty table to collect all data
all_data = table();

% Loop through each CSV file and process the data
for i = 1:length(file_list)
    % Read the CSV file
    file_path = fullfile(folder_path, file_list(i).name);
    file_data = readtable(file_path);

    % Extract channel number from the file name using regular expressions
    tokens = regexp(file_list(i).name, 'Channel(\d+)', 'tokens');
    channel_number = str2double(tokens{1}{1});

    % Initialize an empty table for reshaping data
    reshaped_data = table();

    % Identify unique column prefixes (e.g., "Time_ms", "PosPeak", "NegPeak")
    column_prefixes = {'Time_ms', 'PosPeak', 'NegPeak'};

    % Loop through each prefix and flatten columns
    for j = 1:length(column_prefixes)
        prefix = column_prefixes{j};

        % Find columns with the current prefix
        cols = startsWith(file_data.Properties.VariableNames, prefix);

        % Extract data for these columns and reshape into a single column
        values = table2array(file_data(:, cols));
        values = values(:); % Flatten into a column

        % Remove NaN values that result from reshaping
        values = values(~isnan(values));

        % Create a new column in the reshaped table for each prefix
        reshaped_data.(prefix) = values;
    end

    % Add channel information to each row of reshaped_data
    reshaped_data.Channel = repmat(channel_number, height(reshaped_data), 1);

    % Append the reshaped data to the combined table
    all_data = [all_data; reshaped_data];
end

% Define the output file name
output_file = fullfile(folder_path, "combined_results_4um_derivate.csv");

% Write the combined data to a new CSV file
writetable(all_data, output_file);
disp(['All results combined and saved to ', output_file]);
