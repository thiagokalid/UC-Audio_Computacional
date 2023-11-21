% Função que realiza a síntese de fala tendo como fonte uma matriz com
% os formantes F1,F2,F3 e F4 e o Pitch F0. As larguras de banda são 
%constantes. O sinal de saída é armazenado num ficheiro .wav.


Fs=11025;
%jmax=max(size(vecFOR));
jmax=70;
janela=256;
f0=100;
DF=1/Fs*janela;			% Duração de uma frame.
Av=ones(1,jmax)*1;

% Attack and decay parameters for the envelope
attack = 2 * round(0.2 * Fs / 256 / 2);
decay = 3 * round(0.3 * Fs / 256 / 2);

% Generate the envelope using a Hanning window
e1 = hanning(attack)';
e2 = hanning(decay)';
envolvente = [e1(1:length(e1) / 2) ones(1, jmax - length(e1) / 2 - length(e2) / 2) e2(length(e2) / 2 + 1:length(e2))];

% Apply the envelope to the amplitude
Av = Av .* envolvente;

% Create a vector of pitch values
jitter_factor = 0.05; % Jitter factor (percentage of variation)
F0=readNumbersFromFile("a-gravado.PitchTier");
F0_jittered = F0 + (jitter_factor * randn(1, jmax));

% formantes da vogal a
F1=ones(1,jmax)*700;
F2=ones(1,jmax)*1200;
F3=ones(1,jmax)*2800;
F4=ones(1,jmax)*3400;

% Shimmer factor (percentage of variation)
shimmer_factor = 0.1;
F1_jittered = F1 + (shimmer_factor * randn(1, jmax));
F2_jittered = F2 + (shimmer_factor * randn(1, jmax));
F3_jittered = F3 + (shimmer_factor * randn(1, jmax));
F4_jittered = F4 + (shimmer_factor * randn(1, jmax));

Ts=1/Fs;
B1=ones(1,jmax)*50*6.25;
B2=ones(1,jmax)*75*6.25;
B3=ones(1,jmax)*100*6.25;
B4=ones(1,jmax)*300*6.25;


e=zeros(1,janela);
aant=400*pi;
bant=5000*pi;
Tant=1/10000;
a=aant*Tant/Ts;
b=bant*Tant/Ts;	% Dependem do falante

db=10;
tau1=1e-3;	% 1 ms.
deltatau=2e-3;	% 2 ms.

saida=zeros(1,(jmax+1)*janela/2);

indmax=round(4500*janela/Fs);
vecSPE=zeros(jmax,indmax);
H=hanning(janela);
%sinal=sinal*Av;
ind=1;
resto=[];
zi=zeros(8,1);
% zi=zeros(2,1);
voz=[];excit=[];

% Define a forma do impulso glotal
aglotal=0.2;			
Bglotal=[0 -aglotal*exp(1)*log(aglotal)];
Aglotal=[1 -2*aglotal aglotal^2];

for j=1:jmax,
 % j
    % Generate the glottal pulse
    [sinalglotal,resto]=fgerimp(Fs,F0_jittered(j),janela,resto);
    % plot(sinal)    
    % Apply the envelope to the glottal pulse
    sinal = sinalglotal * Av(j);
    % Append the glottal pulse to the excitation signal
    excit = [excit sinal];
    % Initialize the filter coefficients
    BBB=[1];
    A=[1];

    for i=1:4,
        if i==1, Form=F1_jittered(j);LB=B1(j); end;
        if i==2, Form=F2_jittered(j);LB=B2(j); end;
        if i==3, Form=F3_jittered(j);LB=B3(j); end;
        if i==4, Form=F4_jittered(j);LB=B4(j); end;
        % Calculate the filter coefficients
        BB=1-2*exp(-LB/2*Ts)*cos(2*pi*Form*Ts)+exp(-LB/2*Ts)^2;
        AA=[1 -2*exp(-LB/2*Ts)*cos(2*pi*Form*Ts) exp(-LB/2*Ts)^2];
        % Update the filter coefficients
        BBB=conv(BB,BBB);
        A=conv(A,AA);
    end;

    [sinal,zf]=filter(BBB,A,sinal,zi);
    zi=zf;
    voz=[voz sinal];
    
end;


voz=voz/abs(max(voz));
voz=detrend(voz);
    
eixtemp=(0:jmax-1)*DF*1000/2;
figure(2);

xx=abs(fft(zeros(1,11025)));
for i=1:1000
    sinalaleatorio=randn(1,11025);
    saidaaleatoria=detrend(filter(BBB,A,sinalaleatorio));
    xx=xx+abs(fft(saidaaleatoria));

end
plot(20*log10((xx(2:round(length(xx)/2)))),'r');grid


%reproduce and save 
soundsc(voz,Fs);
audiowrite('A-created.wav', voz, Fs);