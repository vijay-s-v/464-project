clear;
clc;
% load('Capture-2019.04.02.04.49.47.363.mat')
% load('Capture-2019.04.02.04.49.02.081.mat')
load('Capture-2019.04.02.04.38.27.494.mat');

B = 20e6;
[yf dt] = FilRs(Y,XDelta,B);

%==================  MILESTONE 1  ===========================================

% Capture length in ms
capture_length = size(Y, 1) * XDelta * 1000;
timespan = linspace(0, capture_length, size(Y, 1))';

% Plots for #2
figure(1);
subplot(2, 2, 1);
plot(timespan, real(Y));
title('Real part');
subplot(2, 2, 2);
plot(timespan, imag(Y));
title('Imaginary part');
subplot(2, 2, [3, 4]);
plot(timespan, abs(Y));
title('Absolute value');
xlabel('Time [ms]');

%{
figure(2);
subplot(2, 2, 1);
stem(timespan, real(Y));
title('Real part');
subplot(2, 2, 2);
stem(timespan, imag(Y));
title('Imaginary part');
subplot(2, 2, [3, 4]);
stem(timespan, abs(Y));
title('Absolute value');
xlabel('Time [ms]');
%}

P_avg = sum(abs(Y).^2) / (capture_length / 1000);
P_segment = zeros(1, 100);
segment_size = size(Y, 1) / 100;
for i = 1:100
    segment = Y(i + segment_size * (i - 1) : segment_size * i);
    P_segment(1, i) = sum(abs(segment).^2) / 1e-3;
end

figure(3);
stem(P_segment);
title('Average Power of each Segment');
xlabel('Segment');
ylabel('Average Power');

%==================  MILESTONE 2  ===========================================

Y_fft = fft(Y);
N = length(Y);
df = 1 / (N * XDelta);
max_F = N * df;
% Frequency in MHz
frequency_span = linspace(0, max_F, N)' / 1e6;

figure(4);
subplot(2, 2, 1);
plot(frequency_span, real(Y_fft));
title('Real part');
subplot(2, 2, 2);
plot(frequency_span, imag(Y_fft));
title('Imaginary part');
subplot(2, 2, [3, 4]);
plot(frequency_span, abs(Y_fft));
title('Absolute value');
xlabel('Frequency [MHz]');

% Perform FFT on first two 1ms segments
segment1 = Y(1:N/2);
segment2 = Y(1 + N/2:N);
segment1_fft = fft(segment1);
segment2_fft = fft(segment2);

figure(5);
subplot(2, 2, 1);
plot(timespan(1:N/2), abs(segment1));
title('Segment 1');
xlabel('Time [ms]');
subplot(2, 2, 2);
plot(frequency_span(1:N/2), abs(segment1_fft));
title('FFT of Segment 1');
xlabel('Frequency [MHz]');
subplot(2, 2, 3);
plot(timespan(1+N/2:N), abs(segment2));
title('Segment 2');
xlabel('Time [ms]');
subplot(2, 2, 4);
plot(frequency_span(1:N/2), abs(segment2_fft));
title('FFT of Segment 2');
xlabel('Frequency [MHz]');

%==================  MILESTONE 3  ===========================================

% Extract the channel at 11-24MHz
Y_resample = resample(Y, 3072, 5600);
Y_resample_fft = fft(Y_resample);

first = find(frequency_span, 11e6, 'first');
last = find(frequency_span, 24e6, 'last');
Y_resample_fft = Y_resample_fft(first:last);

%==================  FUNCTIONS  ===========================================

function [y dt] = FilRs(xt,dti,B)
    N = length(xt);
    df = 1/(N*dti);
    Nb = floor(B/(2*df));
    Fil(1:N) = 0;
    Fil(1:Nb)= 1;
    Fil(N-Nb+1:N) = 1;
    Xf = fft(xt);
    Xff = Fil'.*Xf;
    dt = dti*56/30.72;
    y = resample(double(ifft(Xff)),3072,5600);
end