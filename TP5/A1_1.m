% LPC Coefficient Calculation with Sliding Analysis

% Load audio file
filename = 'ASRF20.wav';
[x, fs] = audioread(filename);

% Resample the signal to fs=16kHz
fs_new = 16000;
x_resampled = resample(x, fs_new, fs);

% LPC Analysis Parameters
p = 12; % LPC order
window_size = 30e-3; % Window size in seconds
overlap_ratio = 2/3; % Overlap ratio

% Convert window size and overlap to samples
window_length = round(window_size * fs_new);
overlap_length = round(overlap_ratio * window_length);

% Perform LPC analysis with sliding window
num_samples = length(x_resampled);
num_frames = floor((num_samples - overlap_length) / (window_length - overlap_length));

% Initialize arrays to store LPC coefficients and error signals
lpc_coeffs = zeros(p+1, num_frames);
error_signals = zeros(num_samples, 1);

% Hamming window for LPC analysis
hamming_window = hamming(window_length);

% Perform LPC analysis with sliding window
for i = 1:num_frames
    start_idx = (i-1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;
    
    % Extract the current frame
    x_frame = x_resampled(start_idx:end_idx) .* hamming_window;
    
    % LPC Analysis
    lpc_coeffs(:, i) = lpc(x_frame, p);
    
    % Calculate error signal
    error_signals(start_idx:end_idx) = filter(lpc_coeffs(:, i), 1, x_frame);
end

% Display LPC coefficients and plot the error signal
figure;
subplot(2,1,1);
plot((0:num_frames-1)*(window_length-overlap_length)/fs_new, lpc_coeffs');
title('LPC Coefficients');
xlabel('Time (s)');
ylabel('Coefficient Value');

subplot(2,1,2);
t = (0:length(error_signals)-1)/fs_new;
plot(t, error_signals);
title('Error Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Save LPC coefficients and error signal for later use
save('lpc_analysis_results.mat', 'lpc_coeffs', 'error_signals');
