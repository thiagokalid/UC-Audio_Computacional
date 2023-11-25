%% 
% The main goal of this script is to anayse and compare the acquired signals 
% from MATLAB _Data Acquisition Toolbox_ (daq) and the _Analog Input Recorder 
% (air) App_ 

clc
clear

% Load both signals:
[audio_from_daq, Fs1] = audioread("audio_1.wav");
[audio_from_air, Fs2] = aud ioread("audio_2.wav");
% From the acquired signals compute the time interval:
t_span1 = 0:1/Fs1:((length(audio_from_daq) - 1) / Fs1);
t_span2 = 0:1/Fs2:((length(audio_from_air) - 1) / Fs2);
%% 
% The comparissom will be mainly between through the analysis of the raw signals 
% in time domain and spectrogram. 
% 
% In order to compensate eventual gain differences between the two signals, 
% there will be a multiplicative factor:

% Obtain a gain which translate how one signal has more/less energy than the
% other:
gain = rms(audio_from_daq) ./ rms(audio_from_air);
audio_from_air = audio_from_air * gain;
%% 
% After adjusting the energy, the next step is to plot both signals together.


figure()
plot(t_span1, audio_from_daq, color='blue', LineWidth=.1)
title("Comparisson between acquired signals.")
hold on
plot(t_span2, audio_from_air, color='red', LineWidth=.1)
legend('Data Acquisiton Toolbox', 'Analog Input Recorder')
grid()
xlim([0, t_span1(end)]) 
ylim([-1, 1])
xlabel("Time in seconds")
ylabel("Amplitude in volts")