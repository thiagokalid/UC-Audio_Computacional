% LPC Coefficient Calculation with Sliding Analysis
clear all
clc

% Load audio file
filename = 'ASRF20.wav';
[x, fs] = audioread(filename);

% Resample the signal to fs=16kHz
fs_new = 16000;
x_resampled = resample(x, fs_new, fs);

% LPC Analysis Parameters
window_size = 30e-3; % Window size in seconds
overlap_ratio = 2/3; % Overlap ratio

% Convert window size and overlap to samples
window_length = round(window_size * fs_new);
overlap_length = round(overlap_ratio * window_length);

% Extract the portion of the resampled signal from 2.4 to 2.77 seconds
start_time = 2.3; % seconds
end_time = 3.4; % seconds
t_span = start_time:1/fs_new:end_time;
start_index = round(start_time * fs_new);
end_index = round(end_time * fs_new);
x_to_analyze = x_resampled(start_index:end_index);

% Perform LPC analysis with sliding window
num_samples = length(x_to_analyze);
num_frames = floor((num_samples - overlap_length) / (window_length - overlap_length));

% Hamming window for LPC analysis
hamming_window = hamming(window_length);

% Calculate time vector for the analyzed portion
t_analyzed = (start_index:end_index - 1) / fs_new;

% Pre-emphasis filter
pre_emphasis_coefficient = 0.98;
x_to_analyze_preemph = filter([1, -pre_emphasis_coefficient], 1, x_to_analyze);

% Re-synthesize the signal with the original error signal
reconstructed_signal = zeros(size(x_to_analyze));

% Number of LPC coefficients:
p = 20;

% Initialize arrays to store LPC coefficients and error signals
lpc_coeffs1 = zeros(p+1, num_frames);
lpc_coeffs2 = zeros(p+1, num_frames);
error_signal1 = zeros(num_samples, 1);
error_signal2 = zeros(num_samples, 1);
x_est1 = zeros(num_samples, 1);
x_est2 = zeros(num_samples, 1);

% Perform LPC analysis with sliding window
for j = 1:num_frames
    start_idx = (j-1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;

    % Extract the current frame
    x_frame1 = x_to_analyze_preemph(start_idx:end_idx) .* hamming_window;
    x_frame2 = x_to_analyze(start_idx:end_idx) .* hamming_window;

    % LPC Analysis
    lpc_coeffs1(:, j) = lpc(x_frame1, p);
    lpc_coeffs2(:, j) = lpc(x_frame2, p);

    % Calculate error signal
    error_signal1(start_idx:end_idx) = filter(lpc_coeffs1(:, j), 1, x_frame1);
    error_signal2(start_idx:end_idx) = filter(lpc_coeffs2(:, j), 1, x_frame2);

    % Re-synthesize using the LPC coefficients
    x_est1(start_idx:end_idx) = filter(1, lpc_coeffs1(:, j), error_signal1(start_idx:end_idx));
    x_est2(start_idx:end_idx) = filter(1, lpc_coeffs2(:, j), error_signal2(start_idx:end_idx));

end

% De-emphasis:
x_est_preemph = filter(1, [1, -pre_emphasis_coefficient], x_est1);

figure()
hold on
grid on
plot(t_span, x_to_analyze, color='red')
plot(t_span, x_est_preemph, color='blue')
xlim([start_time, end_time])
xlabel("Time in seconds")
ylabel("Amplitude")
ylim([-1, 1])
legend(["Original", ...
    "Estimated with pre-emphasis."])

audiowrite("error_generated_signal.wav", x_est_preemph, fs_new)


