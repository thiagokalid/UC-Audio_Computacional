clear, clc;

% Read the audio data:
[y_old, Fs_old] = audioread("ASRF24.wav");
% Fs_old is 44.1 kHz
t_old = 0:1/Fs_old:(length(y_old)-1)/Fs_old;
Fs_new = 16000; % Hz


% Resampling process:
y_new = resample(y_old, t_old, Fs_new);
t_new = 0:1/Fs_new:(length(y_new)-1)/Fs_new;

%%

hold on
plot(t_old, y_old, ':o')
plot(t_new, y_new, '*')
legend(["F_{s1}=44.1 kHz", "F_{s2}=16 kHz"])
xlim([0.38, 0.41])
ylim([-1, 1])
grid on
xlabel("Time in seconds")
ylabel("Amplitude in volts")

%% Create the desired window:
% Define the window length and the window slide step in seconds:
window_length_in_seconds = 30e-3;
window_step_in_seconds = 10e-3;

% Compute the equivalent size in number of samples:
window_length_in_samples = window_length_in_seconds * Fs_new;
wstep = window_step_in_seconds * Fs_new;

% Imposes that the window length is even:
if (~mod(window_length_in_samples, 2))
    window_length_in_samples = window_length_in_samples + 1;
end
window = hamming(window_length_in_samples);

%% Define which funcions will be evaluated through the sliding analysis:
zcr = moving_window_analysis(y_new, window, wstep, @(m) zerocrossrate(m));
energy = moving_window_analysis(y_new, window, wstep, @(m) sum(m.^2));

%% Normalize to plot together:
energy = (energy - min(energy))/(max(energy) - min(energy));
zcr = (zcr - min(zcr))/(max(zcr) - min(zcr));

%% Plot the resulted value from the sliding analysis.

figure()
lower_t = 0.3;
upper_t = 1.0;
plot(t_new, y_new);
xlim([lower_t, upper_t]);
ylim([-1, 1])

yyaxis right
hold on
t_step = t_new(1:wstep:end);
plot(t_step, zcr, lineWidth=2, color='r');
plot(t_step, energy, lineWidth=2, color='k');
xlim([lower_t, upper_t]);
ylim([-1, 1])

legend("Signal ampltiude", "Zero-crossing rate", 'Energy')
