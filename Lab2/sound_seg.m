clear, clc;

% Read the audio data:
[y_old, Fs_old] = audioread("ASRF24.wav");

t_old = 0:1/Fs_old:(length(y_old)-1)/Fs_old;
Fs_new = 41000/2; % Hz

%
tmin = 9.6;
tmax = 9.8;
i0_constant = 1e-12;

% Resampling process:
y_new = resample(y_old, t_old, Fs_new);
t_new = 0:1/Fs_new:(length(y_new)-1)/Fs_new;

sinal = y_old;
Fs = Fs_old;

%% Create the desired window:
% Define the window length and the window slide step in seconds:
window_length_in_seconds = 30e-3;
window_step_in_seconds = 10e-3;

% Compute the equivalent size in number of samples:
Fs_new = Fs;
window_length_in_samples = window_length_in_seconds * Fs_new;
wstep = window_step_in_seconds * Fs_new;

% Imposes that the window length is even:
if (~mod(window_length_in_samples, 2))
    window_length_in_samples = window_length_in_samples + 1;
end
window = hamming(window_length_in_samples);
window_rec = rectwin(window_length_in_samples);


%%
tmax = .06;
t = 0:1/Fs:tmax;
y = (.6 + .05 * sin(2*pi*600*t)) .* sin(2*pi*100*t);

figure()
hold on
plot(t * 1e3, y, lineWidth=2)
xlabel("Time in ms")
ylabel("Amplitude")

%
t1 = linspace(-15e-3, 15e-3, length(window));
t2 = linspace(-5e-3, 25e-3, length(window));
t3 = linspace(5e-3, 35e-3, length(window));

plot(t1 * 1e3, window, LineWidth=2, color='k')
plot(t2 * 1e3, window, LineWidth=2, color='b')
plot(t3 * 1e3, window, LineWidth=2, color='r')
legend(["Signal", "Window 1", "Window 2", "Window 3"])
grid
xlim([0, tmax * 1e3])



%%

% A Frequencia de Amostragem deve ser 22050 ou 44100 Hz
if Fs==44100
    r = 2;
    sinal=decimate(sinal,r);
end
Fs = Fs/2;
t_span = 0:1/Fs:(length(sinal)-1)/Fs;

sinal=detrend(sinal);
espacamento=10;
energia=fenerg2(sinal,espacamento); % Energia do sinal
m=fmedia(energia,5,1);energia=fmedia(m,3,1); % Alisamento
LSES=max(energia(1:100));
der=fderivad(sinal);
zero=fzeros2(der,50,espacamento); % Taxa de passagem por zero da derivada
m=fmedia(zero,5,1);zero=fmedia(m,3,1);  % Alisamento
LIZ=min(zero(3:70)); % Limite Inferior da Taxa de Passagem por Zero

car=zeros(round(length(sinal)/(espacamento*10)),2);
s=zeros(round(length(sinal)/(espacamento*10)),1);
decisao=zeros(round(length(sinal)/(espacamento*10)),1); 
t_span_decisao = linspace(0, t_span(end), length(decisao));

for i=0:10:length(zero)-10
  z=0;e=0;
  for j=1:10  
    if zero(i+j)>LIZ, z=z+1; end
    if energia(i+j)<LSES, e=e+1; end
  end
  car(i/10+1,1)=e;   % Energia
  car(i/10+1,2)=z;   % Taxa de Passsagem por zero
end
s=sum(car,2);

% Matriz de decisão com 3 classes
% Decisao= 0->nao definido; 1-> silencio; 3-> excitado por ruido (noise); 4-> vocalizado;

for i=1:length(decisao) 
 if s(i)>=14, decisao(i)=1;
   elseif s(i)<6, decisao(i)=4;
     elseif (car(i,1)<3 && car(i,2)>5), decisao(i)=3;
 end
end
for i=3:length(decisao) 
    if decisao(i-2)==decisao(i) 
        decisao(i-1)=decisao(i); 
    end
end


%%



%
decisao_new = interp1(t_span_decisao, decisao, t_span);

%
y_filt = sinal;
y_filt(decisao_new ~= 4) = 0;



%plot(t_span_decisao, y_decisao .* apply_f0)
tmin = 7.83;
tmax = 8.75;

% figure()
% subplot(2, 1, 1)
% 
% 
% 
% plot(t_span, sinal, 'k')
% xlim([tmin, tmax])
% ylim([-.73, .73])
% title("Original signal")
% xlabel("Time in s")
% ylabel("Amplitude")
% grid
% 
% subplot(2, 1, 2)


figure()
x0=10;
y0=10;
width=600;
height=250;
set(gcf,'position',[x0,y0,width,height])
plot(t_span, y_filt, 'k')
% title("Signal containing only voiced parts and computed fundamental frequency (f0).")
ylim([-1, 1])
ylabel("Amplitude")
hold on
yyaxis right

f0_min = 70;
f0_max = 250;
f0 = moving_window_analysis(y_filt, t_span_decisao, window, wstep, @(t, x) find_f0(x, Fs, f0_min, f0_max));
t_span_f0 = t_span(1:wstep:end);

f0(decisao_new(1:wstep:end) ~= 4) = 0;


plot(t_span_f0, f0, '*', MarkerSize=4, color='#D95319')
ylabel("f0 in Hz")
xlabel("Time in s")
grid
ylim([50, 300]);
xlim([tmin, tmax])

%%

[output, y_cropped] = play_segment(t_span, sinal, tmin, tmax, Fs);
% output.play()
% audiowrite("ASRF24_estrangeiros.wav", y_cropped, Fs)