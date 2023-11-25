clear, clc, close all;
Fs = 16000;

tmax = 1;
t = 0:1/Fs:tmax;
f = 1000 * t;
y = sin(2*pi.* f .* t);

%% Create the desired window:
% Define the window length and the window slide step in seconds:
window_length_in_seconds = 30e-3;
window_step_in_seconds = 10e-3;

% Compute the equivalent size in number of samples:
window_length_in_samples = window_length_in_seconds * Fs;
wstep = window_step_in_seconds * Fs;

% Imposes that the window length is even:
if (~mod(window_length_in_samples, 2))
    window_length_in_samples = window_length_in_samples + 1;
end
window = hamming(window_length_in_samples);
window_rec = rectwin(window_length_in_samples);
%% Define which funcions will be evaluated through the sliding analysis:
zcr = moving_window_analysis(y, t, window, wstep, @(t, x) zerocrossrate(x));


%% Plot the results:

plot(t * 1e3, y)
ylim([-2, 2])
ylabel("Amplitude")

yyaxis right
plot(t(1:wstep:end) * 1e3, zcr)
ylabel("Zero-crossing rate")
xlabel("Time in ms")
xlim([0, .3 * 1e3])
legend(["y(t)=sin(2\pi f(t) t) where f(t)=1000t", "Zero-crossing rate"])
