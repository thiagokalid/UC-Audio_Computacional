% Re-synthesize the signal with the original error signal
reconstructed_signal = zeros(size(x_to_analyze));

for j = 1:num_frames
    start_idx = (j-1) * (window_length - overlap_length) + 1;
    end_idx = start_idx + window_length - 1;

    % Extract the current frame
    x_frame = x_to_analyze_preemph(start_idx:end_idx) .* hamming_window;

    % Re-synthesize using the LPC coefficients and error signal
    reconstructed_frame = filter(1, lpc_coeffs(:, j), error_signals(start_idx:end_idx));
    
    % Overlap and add
    reconstructed_signal(start_idx:end_idx) = reconstructed_signal(start_idx:end_idx) + reconstructed_frame;
end

% De-emphasis at the end
reconstructed_signal = filter(1, [1, -pre_emphasis_coefficient], reconstructed_signal);

% Plot the original and re-synthesized waveforms
figure;
subplot(2, 1, 1);
plot(t_analyzed, x_to_analyze, 'b', 'LineWidth', 1);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t_analyzed, reconstructed_signal, 'r', 'LineWidth', 1);
title('Re-synthesized Signal');
xlabel('Time (s)');
ylabel('Amplitude');
