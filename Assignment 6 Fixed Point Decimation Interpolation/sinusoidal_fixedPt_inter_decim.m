%% Input signal for Decimation-Interpolation

n = 0:1:95;
f0 = 100;
f1 = 200;
f2 = 300;
fs = 2400;

xn = sin(2*pi*f0*n/fs) + 0.5*sin(2*pi*f1*n/fs) + 0.6*sin(2*pi*f2*n/fs);

%% Performing floating point operations

M = 2;
L = 2;

fc = 600;
r = 2400;

deciOut = decimator(M, fc, r, 101, xn);
yOut = interpolater(L, fc, r, 101, deciOut);

%% Performing fixed point operations

% format : Q(16,11) => W = 16bits | Q = 11 bits | I = 4 bits | Sign Bits = 1
Q_vals = 1:16;
errors = zeros(length(Q_vals), 1);

for i = 1:length(Q_vals)
    Q = Q_vals(i);
    deciOut_fixed = decimator_fixed(M, fc, r, 101, xn, Q);
    interOut_fixed = interpolater_fixed(L, fc, r, 101, deciOut_fixed, Q);
    yOut_fixed = convert_to_float(interOut_fixed, Q);
    err = calError(yOut, yOut_fixed);
    errors(i) = err;
end

plot(Q_vals, errors);
title('Error vs Q Format')
figure;

[min_error, opt_Q] = min(errors);

% plotting for minimum error Q 

deciOut_fixed = decimator_fixed(M, fc, r, 101, xn, opt_Q);
interOut_fixed = interpolater_fixed(L, fc, r, 101, deciOut_fixed, opt_Q);
yOut_fixed = convert_to_float(interOut_fixed, opt_Q);

[rms_err, errs] = calError(xn, yOut);
fprintf("Decimation-Interpolation Error for M=L=2 : %.4f \n", rms_err);

subplot(3,1,1);
stem(n, xn);
title('Input Waveform')
subplot(3,1,2);
stem(n, yOut_fixed);
title('Output Waveform')
subplot(3,1,3);
stem(n, errs);
title('Errors (M=L=2)')

%% M=L=4

M = 4;
L = 4;
fc = 300;
r = 2400;

deciOut = decimator(M, fc, r, 101, xn);
yOut = interpolater(L, fc, r, 101, deciOut);

% plotting for minimum error Q 

deciOut_fixed = decimator_fixed(M, fc, r, 101, xn, opt_Q);
interOut_fixed = interpolater_fixed(L, fc, r, 101, deciOut_fixed, opt_Q);
yOut_fixed = convert_to_float(interOut_fixed, opt_Q);

[rms_err, errs] = calError(xn, yOut);
fprintf("Decimation-Interpolation Error for M=L=4 : %.4f \n", rms_err);

subplot(3,1,1);
stem(n, xn);
title('Input Waveform')
subplot(3,1,2);
stem(n, yOut_fixed);
title('Output Waveform')
subplot(3,1,3);
stem(n, errs);
title('Errors (M=L=4)')

%% M=L=8

M = 8;
L = 8;
fc = 150;
r = 2400;

deciOut = decimator(M, fc, r, 101, xn);
yOut = interpolater(L, fc, r, 101, deciOut);

% plotting for minimum error Q 

deciOut_fixed = decimator_fixed(M, fc, r, 101, xn, opt_Q);
interOut_fixed = interpolater_fixed(L, fc, r, 101, deciOut_fixed, opt_Q);
yOut_fixed = convert_to_float(interOut_fixed, opt_Q);

[rms_err, errs] = calError(xn, yOut);
fprintf("Decimation-Interpolation Error for M=L=4 : %.4f \n", rms_err);

subplot(3,1,1);
stem(n, xn);
title('Input Waveform')
subplot(3,1,2);
stem(n, yOut_fixed);
title('Output Waveform')
subplot(3,1,3);
stem(n, errs);
title('Errors (M=L=8)')
%% functions definations

function [rms_error, errs] = calError(yFloat, yFixed)
    errs = yFloat - yFixed;
    rms_error = mean(errs.^2);
end

%conversion & operation functions

function arr = convert_to_fixed(float_arr, Q)
    arr = float_arr .* power(2,Q);
    arr = int16(arr);
end

function arr = convert_to_float(fixed_arr, Q)
    arr = double(fixed_arr) ./ power(2,Q);
end

function arr = fixedConv(x, h, Q)
    arr = conv(x, h);
    arr = convert_to_fixed(arr./power(2,2*Q), Q);
end

% fixed point functions

function y = interpolater_fixed(L, fc, fs, N, x, Q)
    xu = upsample(x,L);
    xuLen = length(xu);
    h = convert_to_fixed((L.*finalFilter(fc, fs, N)), Q);
    hlen = length(h);
    y = fixedConv(xu, h, Q);
    y = y((hlen-1)/2+1 : (hlen-1)/2+xuLen);
end

function y = decimator_fixed(M, fc, fs, N, x, Q)

    % converting values to fixed point
    x = convert_to_fixed(x, Q);
    h = convert_to_fixed(finalFilter(fc, fs, N), Q);

    hlen = length(h);
    xlen = length(x);

    xf = fixedConv(x, h, Q);

    xf = xf((hlen-1)/2+1 : (hlen-1)/2+xlen);
    y = downsample(xf, M);

end

% floating point functions

function y = decimator(M, fc, fs, N, x)
    h = finalFilter(fc, fs, N);
    hlen = length(h);
    xlen = length(x);
    xf = conv(x, h);
    xf = xf((hlen-1)/2+1 : (hlen-1)/2+xlen);
    y = downsample(xf, M);
end

function y = interpolater(L, fc, fs, N, x)
    xu = upsample(x,L);
    xuLen = length(xu);
    h = L.*finalFilter(fc, fs, N);
    hlen = length(h);
    y = conv(xu, h);
    y = y((hlen-1)/2+1 : (hlen-1)/2+xuLen);
end

function h = finalFilter(fc, fs, N)
    w = window(1,N);
    hd = filter(fs, fc, N);
    h = hd.*w;
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