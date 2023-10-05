function [output] = moving_window_analysis(signal, window, slide_step, func)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Extract the lenght of the signal and the window. The windo length is
% expected to be always odd.
signal_length = length(signal);
window_length_in_samples = length(window);
half_window_length = (window_length_in_samples - 1)/2;

% The window will be centered at the first element in the first iteration.
% To achieve this effect, '0's will be added to the signal, and then
% ignored after the result is computed.
padded_signal = zeros(signal_length + 2 * half_window_length, 1);
padded_signal(half_window_length + 1: end - half_window_length) = signal;

% The output has the same number of elements as the original signal.
output_length = length(1:slide_step:signal_length);
output = zeros(output_length, 1);
j = 0;
for i = 1:slide_step:signal_length
    j = j + 1;
    idx_beg = i;
    idx_end = i + 2 * half_window_length;
    windowed_signal = padded_signal(idx_beg:idx_end) .* window;
    output(j) = func(windowed_signal);
   
   % if i == 1
   %     figure()
   %     subplot(2, 1, 1)
   %     hold on
   %     title("Original Signal")
   %     plot(padded_signal(idx_beg:idx_end));
   %     plot(window)
   %      legend(["Signal", "Window"])
   %      subplot(2, 1, 2)
   %      plot(windowed_signal)
   %      title("Windowed Signal")
end
end