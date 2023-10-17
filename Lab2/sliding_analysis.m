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
window_rec = rectwin(window_length_in_samples);
%% Define which funcions will be evaluated through the sliding analysis:
zcr = moving_window_analysis(y_new, t_new, window, wstep, @(t, x) zerocrossrate(x));
energy = moving_window_analysis(y_new, t_new, window_rec, wstep, @(t, x) trapz(t, abs(x.^2)));
f0_min = 70;
f0_max = 250;
f0 = moving_window_analysis(y_new, t_new, window_rec, wstep, @(t, x) find_f0(x, Fs_new, f0_min, f0_max));

%% Normalize to plot together:
energy_norm = (energy - min(energy))/(max(energy) - min(energy));
zcr_norm = (zcr - min(zcr))/(max(zcr) - min(zcr));

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
% plot(t_step, zcr_norm, lineWidth=2, color='r');
% plot(t_step, energy_norm, lineWidth=2, color='k');
plot(t_step, f0)
title("Window length = " + num2str(window_length_in_seconds*1e3, 2) + "ms" + ...
    "Window step = " + num2str(window_step_in_seconds*1e3, 2) + "ms")
xlim([lower_t, upper_t]);
ylim([-1, 1])

legend("Signal ampltiude", "Zero-crossing rate", 'Energy')

%%
%% Define which funcions will be evaluated through the sliding analysis:
% f0_min = 70;
% f0_max = 250;
% t = 0:1/Fs_new:1;
% x = sin(2*pi*125 * t);
% find_f0(x, Fs_new, f0_min, f0_max)
% f0 = moving_window_analysis(y_new, t, window_rec, wstep, @(t, x) find_f0(x, Fs_new, f0_min, f0_max));

%% Normalize to plot together:
energy_norm = (energy - min(energy))/(max(energy) - min(energy));
zcr_norm = (zcr - min(zcr))/(max(zcr) - min(zcr));
f0_norm = (f0 - min(f0))/(max(f0) - min(f0));

%% Read PRAAT results:
praat_energy = readmatrix('energy.csv');
praat_energy_tspan = linspace(0, 20.848390022675737, length(praat_energy));
praat_f0 = readmatrix('f0.csv');
praat_f0_tspan = linspace(0, 20.848390022675737, length(praat_f0));

%% Plot normalized zero crossing rate, energy and f0
figure()
lower_t = 0.3;
upper_t = 1.0;
plot(t_new, y_new);
xlim([lower_t, upper_t]);
ylim([-1, 1])
plot(t_new, y_new)
yyaxis right
hold on
t_step = t_new(1:wstep:end);
plot(t_step, zcr_norm, lineWidth=2, color='r');
plot(t_step, energy_norm, lineWidth=2, color='k');
plot(t_step, f0_norm, linewidth=2, color='g')
title("Window length = " + num2str(window_length_in_seconds*1e3, 2) + "ms" + ...
    "Window step = " + num2str(window_step_in_seconds*1e3, 2) + "ms")
xlim([0, 5]);
ylim([-1, 1])
legend("Signal ampltiude", "Zero-crossing rate", 'Energy', 'f_0')


%%

figure()
subplot(3,1,1);
plot(t_step, 20*log(energy), color='r')
hold on
grid on
yyaxis right
plot(praat_energy_tspan, praat_energy, color='b')
title("Energyt")
legend(["MATLAB", "PRAAT"])

subplot(3,1,2);
plot(t_step, zcr, color='r')
hold on
grid on
yyaxis right
plot(praat_energy_tspan, praat_energy, color='b')
title("Zero-crossing rate")
legend(["MATLAB", "PRAAT"])

subplot(3,1,3);
plot(t_step, f0, color='r')
hold on
grid on
yyaxis right
plot(praat_energy_tspan, praat_energy, color='b')
title("Fundamental Frequency (f0)")
legend(["MATLAB", "PRAAT"])
