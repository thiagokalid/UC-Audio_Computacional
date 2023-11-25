function [y] = sliding_func(x, window, step, func)
window_length = length(window);
half_window_length = (window_length - 1)/2;
y = zeros(length(x));
x_padded = zeros(length(x) + 2*half_window_length);
x_padded(half_window_length:end - half_window_length) = x;
for i = 1:length(x)
    sliced_x = x_padded(i + half_window_length);
    y(i) = func(sliced_x * window);
end