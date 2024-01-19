clear, close, clc
% Read the audio file
filename = 'Sounds/teste_casa3.wav';
[x, Fs] = audioread(filename);


% Undersample to 16 kHz
targetFs = 16000;
x = resample(x, targetFs, Fs);
Fs = targetFs;
t = 0:1/Fs:(length(x)-1)/Fs;

%% Set dynamic range compressor parameters:
T = -15;
R = 10;
W = 5;
AttackTime = 0.05;
ReleaseTime = 0.4;

maxAttackTime = 200e-3;
maxReleaseTime = 1;

timeConstant = 200e-3; % Time-constant defined according to reference;

MakeUpGain = 0;

%% Apply the dynamic range compressor:
[x_compressed, gain_linear, attackTime, releaseTime] = smart_compressor(x, Fs, T, R, W, maxAttackTime, maxReleaseTime, timeConstant, MakeUpGain);
[x_compressed2, gain_linear2, attackTime2, releaseTime2] = classic_compressor(x, Fs, T, R, W, AttackTime, ReleaseTime, timeConstant, MakeUpGain);

%% Plot the original audio, output of the compressor and gain:

figure()
subplot(2,1,1)
hold on
grid on
plot(t, x);
plot(t, x_compressed)
Tvector = 10^(T/20) + gain_linear .* 0;
plot(t, Tvector, ":k", "LineWidth", 1.5)
plot(t, -Tvector, ":k", "LineWidth", 1.5)
% title("Comparisson between the compressor input and ouput.")
xlim([t(1), t(end)])
ylim([-1.1, 1.1])
legend(["Original", "Compressed", "Threshold"])
ylabel("Signal Amplitude")
xlabel("Time in seconds")

subplot(2,1,2)
plot(t, gain_linear, "Color", "blue", "LineWidth", 2)
% title("Compressor Gain (in linear scale)")
xlim([t(1), t(end)])
ylim([-.1, 1.1])
xlabel("Time in seconds")
grid on

%% Plot the attack time and release time in different axis:

figure()
subplot(2,1,1)
hold on
grid on
plot(t, x);
plot(t, x_compressed)
Tvector = 10^(T/20) + gain_linear .* 0;
plot(t, Tvector, ":k", "LineWidth", 1.5)
plot(t, -Tvector, ":k", "LineWidth", 1.5)
% title("Comparisson between the compressor input and ouput.")
xlim([t(1), t(end)])
ylim([-1.1, 1.1])
legend(["Original", "Compressed", "Threshold"])
ylabel("Signal Amplitude")
xlabel("Time in seconds")

subplot(2,1,2)
hold on
grid on
plot(t, attackTime * 1e3, 'r', "LineWidth", 1.5)
plot(t, releaseTime * 1e3, 'b', "LineWidth", 1.5)
legend(["Attack Time", "Release Time"])

ylabel("Time in ms")
xlabel("Time in s")
xlim([t(1), t(end)])

%%
figure()
subplot(2,1,1)
plot(t, x, "color", "blue")
grid on
hold on
plot(t, x_compressed2, "LineWidth", 1.5, "color", "red")
plot(t, x_compressed, "LineWidth", 1.5, "color", "magenta")
plot(t, Tvector, ":k", "LineWidth", 1.5)
plot(t, -Tvector, ":k", "LineWidth", 1.5)
legend(["Original", "Classic", "Smart"])
xlim([t(1), t(end)])

subplot(2,1,2)
plot(t, gain_linear, "LineWidth", 1.5, "color", "magenta")
grid on
hold on
plot(t, gain_linear2, "LineWidth", 1.5, "Color", "red")
legend(["Smart", "Classic"])
xlim([t(1), t(end)])

