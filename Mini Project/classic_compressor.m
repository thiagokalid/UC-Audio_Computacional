function [audioOut, gainOut, AttackTime, ReleaseTime] = classic_compressor(audioIn, Fs, threshold, ratio, kneeWidth, AttackTime, ReleaseTime, timeConstant, MakeUpGain)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x = audioIn;
T = threshold;
R = ratio;
W = kneeWidth;

% Transform to log-scale (dB):
x_db = 20 * log10(abs(x));

% Initiallise vectors that will be used:
x_sc = zeros(size(x_db)); 
gm = zeros(size(x_db)); 
g = ones(size(x_db)); 
y = zeros(size(x_db)); 
gs = zeros(size(x_db));

for n = 2:length(x_db)
    % Gain computing:
    if x_db(n) < (T - W/2)
        x_sc(n) = x_db(n);
    elseif (T - W/2) < x_db(n) && x_db(n) < (T + W/2)
        x_sc(n) = x_db(n) + ((1/R - 1)*(x_db(n) - T + W/2)^2)/(2*W);
    else
        x_sc(n) = T + (x_db(n) - T)/R;
    end

    gc(n) = x_sc(n) - x_db(n);

    % Smoothen the value of the gain according to attackTime and releaseTime:
    alpha_a = exp(-log(9)/(AttackTime * Fs)); % Constant
    alpha_r = exp(-log(9)/(ReleaseTime * Fs)); % Constant

    if gc(n) <= gs(n-1)
        gs(n) = alpha_a * gs(n-1) + (1-alpha_a) * gc(n);
    else
        gs(n) = alpha_r * gs(n-1) + (1-alpha_r) * gc(n);
    end

    % Compute the make-up gain.
    M = 0;

    % Apply make-up gain:
    gm(n) = gs(n) + M;

    % Transform back the gain to linear scale:
    g(n) = 10^(gm(n)/20);
    
    % Apply the gain:
    y(n) = x(n) * g(n);

end
audioOut= y;
gainOut = g;
end