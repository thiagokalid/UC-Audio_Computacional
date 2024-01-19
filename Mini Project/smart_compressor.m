function [audioOut, gainOut, AttackTime, ReleaseTime] = smart_compressor(audioIn, Fs, threshold, ratio, kneeWidth, maxAttackTime, maxReleaseTime, timeConstant, MakeUpGain)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x = audioIn;
T = threshold;
R = ratio;
W = kneeWidth;

% Transform to log-scale (dB):
x_db = 20 * log10(abs(x));

% Initialise vectors the cases where the compressing is supposed to work:
x_sc = zeros(size(x_db)); 
gs = zeros(size(x_db));
gc = zeros(size(x_db));
y_peak = zeros(size(x_db));
y_rms = zeros(size(x_db));
y_c = zeros(size(x_db));
c_dev = zeros(size(x_db));
c_makeup = zeros(size(x_db));
c = zeros(size(x_db));
y = zeros(size(x_db));
g = ones(size(x_db));
AttackTime = zeros(size(x_db));
ReleaseTime = zeros(size(x_db));

% Attack-time and Release-time related constants:
alpha = exp(-1/(timeConstant * Fs));
c_est = 10.^(T/20)*(1-1/10.^(R/20))/2;

for n = 2:length(x_db)
    if x_db(n) < (T - W/2)
        x_sc(n) = x_db(n);
    elseif (T - W/2) < x_db(n) && x_db(n) < (T + W/2)
        x_sc(n) = x_db(n) + ((1/R - 1)*(x_db(n) - T + W/2)^2)/(2*W);
    else
        x_sc(n) = T + (x_db(n) - T)/R;
    end

    % Smoothen the value of the gain according to attackTime and releaseTime:
    gc(n) = x_sc(n) - x_db(n);
    
    % Peak detector:
    y_peak(n) = sqrt(max([x(n)^2, ...
                  alpha * y_peak(n-1)^2 + (1 - alpha) * abs(x(n)^2)...
                  ]));

    % RMS detector:
    y_rms(n) = sqrt(alpha * y_rms(n - 1)^2 + (1-alpha) * x(n)^2);

    % Compute the crest-factor:
    y_c(n) = y_peak(n)/y_rms(n);
    
    % Compute the attack-time and release-time:
    AttackTime(n) = 2 * maxAttackTime ./ (y_c(n).^2);
    ReleaseTime(n) = 2 * maxReleaseTime ./ (y_c(n).^2) - AttackTime(n);
    
    % Make sure that they are at most equal to their maximum user-defined
    % value:
    if AttackTime(n) > maxAttackTime || isnan(AttackTime(n))
        AttackTime(n) = maxAttackTime;
    end
    if ReleaseTime(n) > maxReleaseTime || isnan(ReleaseTime(n))
        ReleaseTime(n) = maxReleaseTime;
    end
    
    % Smoothen the value of the gain according to the instantaneous 
    % attackTime and releaseTime:
    alpha_a = exp(-log(9)/(AttackTime(n) * Fs));
    alpha_r = exp(-log(9)/(ReleaseTime(n) * Fs));

    if gc(n) <= gs(n-1)
        gs(n) = alpha_a * gs(n-1) + (1-alpha_a) * gc(n);
    else
        gs(n) = alpha_r * gs(n-1) + (1-alpha_r) * gc(n);
    end
    c(n) = 10.^(gs(n)/20);
    c_dev(n) = alpha * c_dev(n-1) + (1- alpha) * (gs(n) - c_est);
    c_makeup(n) = -(c_dev(n) + c_est);
    
    % Compute the make-up gain.
    % M = -(c_dev + c_est);
    M = 0; % We decided to use 0 dB for make-up gain just to simplify the
    % later comparissons between methods.
    gm(n) = gs(n) + M;

    % Transform back the gain to linear scale:
    g(n) = 10^(gm(n)/20);

    % Apply the gain:
    y(n) = x(n) * g(n);    
end
audioOut = y;
gainOut = g;
end