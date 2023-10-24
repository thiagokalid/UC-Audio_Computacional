clear, clc, close all;

% Read the audio data:
[y_old, Fs_old] = audioread("ASRF24.wav");
% Fs_old is 44.1 kHz
t_old = 0:1/Fs_old:(length(y_old)-1)/Fs_old;
Fs_new = 16000; % Hz

%
tmin = 10.44;
tmax = 16.12;
i0_constant = 1e-12;

% Resampling process:
y_new = resample(y_old, t_old, Fs_new);
t_new = 0:1/Fs_new:(length(y_new)-1)/Fs_new;
% plot(t_new, y)


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
energy = moving_window_analysis(y_new, t_new, window_rec, wstep, @(t, x) trapz(t, x.^2));
energy_norm = (energy - min(energy))/(max(energy) - min(energy));
zcr_norm = (zcr - min(zcr))/(max(zcr) - min(zcr));

%% Estimate signal/silence regions:

noise_energy = mean(energy(1:35));
energy_threshold = (max(energy) - noise_energy) * .001;
silence = energy < energy_threshold;

%% Plot the results:

figure()
% subplot(2,1,1)
plot(t_new, y_new)
ylabel("Amplitude")
yyaxis right
plot(t_new(1:wstep:end), silence, 's', MarkerSize=2)
xlim([tmin, tmax])
ylabel("Classification label. 0=Signal, 1=Silence")
ylim([-1.1, 1.1])
grid on

%% Play the audio clip and save it in .wav format:

[output, y_cropped] = play_segment(t_new, y_new, tmin, tmax, Fs_new);
output.play()
audiowrite("ASRF24_silence.wav", y_cropped, Fs_new)