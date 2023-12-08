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

% Compute the average fundamental frequency from PRAAT:
f0 = mean(readNumbersFromFile("segment_to_analyze.PitchTier"));

% Input the voiced segments found in PRAAT:
t_beg = [.14798 0.617085] + start_time;
t_end = [.34699 0.986794] + start_time;

% 
f0_vector = zeros(size(t_analyzed));
f0_vector(t_beg(1) < t_analyzed & t_analyzed < t_end(1)) = f0;
f0_vector(t_beg(2) < t_analyzed & t_analyzed < t_end(2)) = f0;


% Define the number of LPC coefficients:
p = 20;

% Initialize arrays to store LPC coefficients and error signals
lpc_coeffs = zeros(p+1, num_frames);
error_signal = zeros(num_samples, 1);

% Perform LPC analysis with sliding window
for j = 1:num_frames
    start_idx = (j-1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;

    % Extract the current frame
    x_frame = x_to_analyze_preemph(start_idx:end_idx) .* hamming_window;

    % LPC Analysis
    lpc_coeffs(:, j) = lpc(x_frame, p);

    % Calculate error signal
    error_signal(start_idx:end_idx) = filter(lpc_coeffs(:, j), 1, x_frame);

    % Calculate the variance of the error:
    signal_power = sum(x_frame.^2);
    
    segmentSize = end_idx - start_idx + 1;
    if (f0_vector(start_idx) > 0)
        % If it is a voiced segment:
        step = round(fs_new/f0_vector(start_idx));
        src_signal = zeros(segmentSize, 1);
        idx = offset+1:step:segmentSize;
        offset = step + idx(end) - segmentSize;
        src_signal(idx) = sqrt(step);
    else
        % If it is a unvoiced segment:
        offset = 0;
        src_signal = randn(segmentSize, 1);

    end

    % Re-synthesize using the LPC coefficients and noise/pulse signal:
    reconstructed_signal(start_idx:end_idx) = reconstructed_signal(start_idx:end_idx) + ...
        filter(1, [1; lpc_coeffs(2:end, j)], signal_power * src_signal);

end


reconstructed_signal = filter(1, [1, -pre_emphasis_coefficient], reconstructed_signal);

soundsc(reconstructed_signal, fs_new)
num2str(sum((x_to_analyze - reconstructed_signal).^2), 4)

figure()
hold on
grid on
plot(t_span, reconstructed_signal/max(reconstructed_signal), color='blue')
plot(t_span, x_to_analyze, color='red')
plot(t_analyzed, f0_vector/max(f0_vector))
xlim([start_time, end_time])
ylim([-1, 1])
xlabel("Time in seconds")
ylabel("Amplitude")
legend(["Reconstructed", "Original", "1=Voiced, 0=Unvoiced)"])

audiowrite("generated_excitation_signal.wav", reconstructed_signal, fs_new)


