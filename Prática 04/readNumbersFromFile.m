function [pitch] = readNumbersFromFile(filename)

% Read numbers from a text file and extract relevant data

% Open the file
file = fopen(filename, 'r');

% Check if the file was opened successfully
if file < 0
    error('Failed to open file');
end

% Initialize empty arrays to store extracted data
line = [];
data = [];

% Read the file line by line
while feof(file) == 0
    line = fgetl(file);
    data = [data; str2double(line)];
end

% Extract start and end times
start_time = data(4);
end_time = data(5);

% Extract odd-numbered pitch values
pitch = data(6:2:end);

% Close the file
fclose(file);