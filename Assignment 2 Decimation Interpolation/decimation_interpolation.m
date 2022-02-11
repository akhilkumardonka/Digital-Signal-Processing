%% Variable Initialization
clc;
close all;

%% Input signal for Decimation-Interpolation

n = 0:1:95;
f0 = 100;
f1 = 200;
f2 = 300;
fs = 2400;

xn = sin(2*pi*f0*n/fs) + 0.5*sin(2*pi*f1*n/fs) + 0.6*sin(2*pi*f2*n/fs);
stem(n, xn);

%% Question 1:
M = 2;
L = 2;
deciOut = decimator(M, 600, 2400, 101, xn);
yOut = interpolater(L, 600, 2400, 101, deciOut);
[err, errors] = calError(xn, yOut);
fprintf("Decimation-Interpolation Error for M=L=2 : %.4f \n", err);

subplot(2,1,1);
stem(n, yOut);
title('Output Waveform')
subplot(2,1,2);
stem(n, errors);
title('Errors')

%% Question 2
M = 4;
L = 4;
deciOut = decimator(M, 300, 2400, 101, xn);
yOut = interpolater(L, 300, 2400, 101, deciOut);
[err, errors] = calError(xn, yOut);
fprintf("Decimation-Interpolation Error for M=L=4 : %.4f\n", err);
subplot(2,1,1);
stem(n, yOut);
title('Output Waveform')
subplot(2,1,2);
stem(n, errors);
title('Errors')

%% Question 3
M = 8;
L = 8;
deciOut = decimator(M, 150, 2400, 101, xn);
yOut = interpolater(L, 150, 2400, 101, deciOut);
[err, errors] = calError(xn, yOut);
fprintf("Decimation-Interpolation Error for M=L=8 : %.4f \n", err)
subplot(2,1,1);
stem(n, yOut);
title('Output Waveform')
subplot(2,1,2);
stem(n, errors);
title('Errors')

%% functions definations

function [error, errs] = calError(x, y)
    errs = y - x;
    error = mean(abs(errs));
end

function y = decimator(M, fc, fs, N, x)
    w = window(1,N);
    hd = filter(fs, fc, N);
    h = hd.*w;
    hlen = length(h);
    xlen = length(x);
    xf = conv(x, h);
    xf = xf((hlen-1)/2+1 : (hlen-1)/2+xlen);
    y = downsample(xf, M);
end

function y = interpolater(L, fc, fs, N, x)
    xu = upsample(x,L);
    xuLen = length(xu);
    w = window(1,N);
    hd = filter(fs, fc, N);
    h = L.*(hd.*w);
    hlen = length(h);
    y = conv(xu, h);
    y = y((hlen-1)/2+1 : (hlen-1)/2+xuLen);
end

function hd = filter(fs, fc, N)
    hd = zeros(1, N); % preallocation
    range = -(N-1)/2:1:(N-1)/2;
    wc = (2*pi*fc)/fs; % digital frequency
    for n = 1:length(range)
        if range(n)~=0
            hd(n) = sin(wc*range(n))/(pi*range(n));
        else
            hd(n) = wc/pi;
        end
    end
end

function win = window(style, N)
    winRange = 0:1:N-1;
    if style == 1
        win = 0.54 - 0.46*cos(2*pi.*winRange/(N-1)); %hamming window
    elseif style == 2
        win = 0.5 - 0.5*cos(2*pi.*winRange/(N-1)); %hanning window
    elseif style == 3
        win = 0.42 - 0.5*cos(2*pi.*winRange/(N-1)) + 0.08*cos(4*pi.*winRange/(N-1)); %blackman window
    end
end

