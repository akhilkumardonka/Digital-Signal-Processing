clc;
close all;
%% getting audio samples and plotting spectrogram

[d,r] = audioread('msmn1.wav');

%% Decimation and Interpolation | Modify M, L 

M = 2;
L = 2;
fc = r / (2*M);
deciOut = decimator(M, fc, r, 101, d);
% audiowrite('decimated_output.wav',deciOut, round(r/M));
df = interpolater(L, fc, r, 101, deciOut);
soundsc(df,r);
% audiowrite('interpolated_output.wav',df, r);

%% Analyzing output

subplot(3,1,1);
specgram(d,1024,r);
title('Original Spectrogram')
subplot(3,1,2);
specgram(deciOut,1024,r/M);
title('Decimated Spectrogram')
subplot(3,1,3);
specgram(df,1024,r);
title('Decimated Interpolated Spectrogram')

%% Functions in Use

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

