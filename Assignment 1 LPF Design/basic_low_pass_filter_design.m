%% Variable initialization and filter parameters
clc;
close all;

wc = pi/2;
fc = 400;
N=39;
fs = (2*pi*fc)/wc;

%% array initialization
range = -(N-1)/2:1:(N-1)/2;

%% filter designing
hd = zeros(1,N);

for n = 1:length(range)
    if range(n)~=0
        hd(n) = sin(wc*range(n))/(pi*range(n));
    else
        hd(n) = wc/pi;
    end
end

%% plotting spectrum(fft)
fvtool(hd);

%% Windowing (hamming)
winRange = 0:1:N-1;
w = 0.54 - 0.46*cos(2*pi.*winRange/(N-1)); %hamming window
h =  hd.*w; % applying window
fvtool(h);