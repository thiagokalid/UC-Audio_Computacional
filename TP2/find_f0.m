function [f0] = find_f0(signal, Fs, f0_min, f0_max)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[autocrr, lags] = xcorr(signal);

% Compute the search window:
N = length(signal);
lagmin = ceil(N + Fs /f0_max);
lagmax = ceil(N + Fs / f0_min);

% Apply the window to the autocorrelation function:
windowed_autocrr = autocrr(lagmin:lagmax);
windowed_lags = lags(lagmin:lagmax);

%stem(windowed_lags, windowed_autocrr)

% Find the peak which is the fundamental period:
[peak_value, peak_position] = max(windowed_autocrr);
shift_value = windowed_lags(peak_position);
fundamental_period = shift_value / Fs;
f0 = 1/fundamental_period;

end

% fs = 16100;
% t = 0:1/fs:1;
% f0_min = 70;
% f0_max = 250;
% 
% x1 = sin(2*pi*125 * t);
% 
% N = length(x1);
% lagmin = ceil(N + fs /f0_max);
% lagmax = ceil(N + fs / f0_min);
% [autocrr, lags] = xcorr(x1);
% windowed_autocrr = autocrr(lagmin:lagmax);
% windowed_lags = lags(lagmin:lagmax);
% stem(windowed_lags, windowed_autocrr)
% [peak_value, peak_position] = max(windowed_autocrr)
% shift_value = windowed_lags(peak_position)
% fundamental_period = shift_value / fs
% f0 = 1/fundamental_period;
% disp(f0)