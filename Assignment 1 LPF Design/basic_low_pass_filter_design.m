%% Variable initialization and filter parameters
clc;
close all;

%% filter designing | filter(fs, fc, N)
hd = filter(1600, 400, 39);
fvtool(hd);

%% Windowing (hamming)

% window(style, N) | style : 1-Hamming, 2-Hanning, 3-Blackman
w = window(3,39);
h =  hd.*w; % applying window
fvtool(h);

%% Functions

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