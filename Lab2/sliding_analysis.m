clear, clc;

% Read the audio data:
[y_old, Fs_old] = audioread("ASRF24.wav");
% Fs_old is 44.1 kHz
t_old = 0:1/Fs_old:(length(y_old)-1)/Fs_old;
Fs_new = 16000; % Hz

%
tmin = 7.83;
tmax = 8.75;
i0_constant = 1e-12;

% Resampling process:
y_new = resample(y_old, t_old, Fs_new);
t_new = 0:1/Fs_new:(length(y_new)-1)/Fs_new;

%% Plot a clip of the downsampled signal:

hold on
plot(t_old, y_old, ':o')
plot(t_new, y_new, '*')
legend(["F_{s1}=44.1 kHz", "F_{s2}=16 kHz"])
xlim([0.38, 0.41])
ylim([-1, 1])
grid on
xlabel("Time in seconds")
ylabel("Amplitude in volts")
title("Downsampling of the original siginal.")

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
energy = moving_window_analysis(y_new, t_new, window, wstep, @(t, x) trapz(t, x.^2));
f0_min = 70;
f0_max = 250;
f0 = moving_window_analysis(y_new, t_new, window, wstep, @(t, x) find_f0(x, Fs_new, f0_min, f0_max));

%% Normalize to plot together:
energy_norm = (energy - min(energy))/(max(energy) - min(energy));
zcr_norm = (zcr - min(zcr))/(max(zcr) - min(zcr));
f0_norm = (f0 - min(f0))/(max(f0) - min(f0));

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
plot(t_step, energy_norm, LineWidth=2)
plot(t_step, zcr_norm, LineWidth=2)
plot(t_step, f0_norm, 's')
title("Signal and its short-term metrics. Window length = " + num2str(window_length_in_seconds*1e3, 2) + "ms " + ...
    "Window step = " + num2str(window_step_in_seconds*1e3, 2) + "ms")
xlim([lower_t, upper_t]);
ylim([-1, 1])
ylabel("Normalized value")

legend("Signal ampltiude", 'Energy (Normalized)', "Zero-crossing rate (Normalized)", 'f0 (Normalized)')




%% Read PRAAT results:
praat_energy = readmatrix('energy.csv');
praat_energy_tspan = linspace(0, 20.848390022675737, length(praat_energy));
praat_f0 = readmatrix('f0.csv');
praat_f0_tspan = linspace(0, 20.848390022675737, length(praat_f0));

%% Plot normalized zero crossing rate, energy and f0
figure()
x0=10;
y0=10;
width=600;
height=250;
set(gcf,'position',[x0,y0,width,height])

plot(t_new, y_new, 'k');
xlim([tmin, tmax]);
ylim([-1, 1])
ylabel("Amplitude")
xlabel("Time in s")

yyaxis right
hold on
t_step = t_new(1:wstep:end);
plot(t_step, f0, '*', MarkerSize=4, color='#D95319')
title("Signal and its fundamental frequncy (f0).")
ylabel("Fundamental Frequency in Hz")
grid

%%

figure()
x0=10;
y0=10;
width=600;
height=250;
set(gcf,'position',[x0,y0,width,height])
plot(praat_energy_tspan, praat_energy, color='b')
ylabel("PRAAT Intensity in dB");
hold on
grid on

yyaxis right
ylabel("MATLAB Intensity in dB")
xlabel("Time in s")
plot(t_step, 20*log(energy), color='r')
legend(["MATLAB", "PRAAT"])
xlim([0, 20])
title("Comparisson between MATLAB and PRAAT computed short-term energy")

%%
[output, y_cropped] = play_segment(t_new, y_new, tmin, tmax, Fs_new);
output.play()
% audiowrite("ASRF24_5s.wav", y_cropped, Fs)