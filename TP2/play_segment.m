function [audioplayer_obj, signal_cropped] = play_segment(t, y, tmin, tmax, Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[m, idx_beg] = min((t - tmin).^2);
[m, idx_end] = min((t - tmax).^2);
signal_cropped = y(idx_beg:idx_end);
audioplayer_obj = audioplayer(signal_cropped, Fs);
end