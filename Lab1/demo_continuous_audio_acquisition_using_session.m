% Firstly, we clear the command window and all old variables:
clc
clear

%% Firstly we will be configuring properly the acquisition system.
% On the given case, the acquisiton system is audio input and the API
% used to control this system is the 'directsound'.

% Check all available devices:
d = daq.getDevices

% Chooses the system default audio input device:
dev = d(1)

%% Create an audio session
% Create a session with |directsound| as the vendor and add an audio input
% channel to it. The purpose of the sesison object is to store data and
% many relavant metadata such as sampling frequency.

s = daq.createSession('directsound');
ch1 =   addAudioInputChannel(s, dev.ID, 1:1);

% Set the acquisiton duration to be:
s.DurationInSeconds = 10 % in seconds
% Obtain the sampling frequency. It importance in this script is in the
% saving of the acquired data in .wave format:
Fs = s.Rate

%% 
% Prepare session for continuous operation. The parameter "isContinous"
% sets the necessity or not of calling a "stop(s)" after the "start(s)".
% If this property is false, the operation will be held after the time of
% s.DurationInSeconds.
s.IsContinuous = false

%% 
% Set up the plot for an FFT of the live input.

hf = figure(1);  
hp = plot(zeros(1000,1));  
T = title('Discrete FFT Plot');
xlabel('Frequency (dB Hz)')
ylabel('|Y(f)|dB')
grid on;

%% Add |DataAvailable| listener
% Listener updates the figure with the FFT of the live input signal.

plotFFT = @(src, event) helper_continuous_fft(event.Data, src.Rate, hp, false);
hl = addlistener(s, 'DataAvailable', plotFFT);

%% Start acquisition
% The function startForeground not only starts the acquisiton, but holds the
% run in this line during the acquisition time span.
[y, x, trigger_time] = startForeground(s);

%%
delete(hl);

%% Plot the acquired data and its spectrogram:
figure(2)
subplot(2,1,1)
plot(x, y)
xlim([0, s.DurationInSeconds])
ylim([-1.05, 1.05])
title("Acquired signal in time domain.")
xlabel("Time /[s]")
ylabel("Amplitude /[V]")
grid

subplot(2,1,2)

window_time_span = .06;
L = window_time_span * Fs;
window = hamming(L);
colormap gray;
spectrogram(y, window, "yaxis",[],[],Fs);
title("Spectogram")
xlabel("Time /[s]")
ylabel("Frequency /[kHz]")
ylim([0, 1])

%% Save the output
% Save the vector containing the audio data in a .wave file.
audiowrite("audio_1.wav", y, Fs)