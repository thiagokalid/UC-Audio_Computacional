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

% Explore the influence of p on the characteristics of the error signal
p_values = 0:1:20; % Experiment with different LPC orders
num_p_values = length(p_values);
rmse_values = zeros(num_p_values, 1);
residual_energy = zeros(num_p_values, 1);

% Pre-emphasis filter
pre_emphasis_coefficient = 0.98;
x_to_analyze_preemph = filter([1, -pre_emphasis_coefficient], 1, x_to_analyze);

% Re-synthesize the signal with the original error signal
reconstructed_signal = zeros(size(x_to_analyze));

for i = 1:num_p_values
    p = p_values(i);

    % Initialize arrays to store LPC coefficients and error signals
    lpc_coeffs = zeros(p+1, num_frames);
    error_signals = zeros(num_samples, 1);

    % Perform LPC analysis with sliding window
    for j = 1:num_frames
        start_idx = (j-1) * (window_length - overlap_length) + 1;
        end_idx = start_idx + window_length - 1;

        % Extract the current frame
        x_frame = x_to_analyze_preemph(start_idx:end_idx) .* hamming_window;

        % LPC Analysis
        %lpc_coeffs(:, j) = lpc(x_frame, p);
        lpc_coeffs(:, j) = lpc(x_frame, p);

        % Calculate error signal
        error_signals(start_idx:end_idx) = filter(lpc_coeffs(:, j), 1, x_frame);

        % Re-synthesize using the LPC coefficients and error signal
        reconstructed_frame = filter(1, lpc_coeffs(:, j), pulse_train);

        % Overlap and add
        reconstructed_signal(start_idx:end_idx) = reconstructed_signal(start_idx:end_idx) + reconstructed_frame;
    end
    % Calculate RMSE
    rmse_values(i) = sqrt(mean((error_signals - filter(lpc_coeffs(:, j), 1, x_to_analyze)).^2));

    % Calculate Residual Energy
    residual_energy(i) = sum(error_signals.^2);


end

reconstructed_signal = filter(1, [1, pre_emphasis_coefficient], x_to_analyze);

sound(x_to_analyze, fs_new)
sleep(1)
sound(reconstructed_signal, fs_new)

figure()
hold on
grid on
plot(t_span, reconstructed_signal/max(reconstructed_signal))
plot(t_span, x_to_analyze)
xlim([start_time, end_time])
ylim([-1, 1])
xlabel("Time in s")
ylabel("Amplitude")
legend(["Reconstruido", "Original"])


