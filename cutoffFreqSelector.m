% x = dffAllConds{4, 1}(:, 10);
x = dFFThisHolo;
% x = meanFiltDffAllConds5:, 1);

% Sampling frequency
fs = imagingFreq;

% Number of samples
N = length(x);

% Perform FFT
X = fft(x);

% Frequency vector (up to Nyquist)
f = (0:N-1) * (fs/N);

% Normalize and compute magnitude
X_mag = abs(X)/N;

% Plot only the positive frequencies (up to Nyquist)
half_N = floor(N/2);
figure(101);
plot(f(1:half_N), X_mag(1:half_N));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of the Signal');
grid on;

figure(1);
clf
cutOffFreq = 20;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt1 = x;
filtdff = filter(b, a, filt1);
plot(filtdff);
% ylim([-0.1 0.2])

figure(2);
clf
cutOffFreq = 30;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt2 = x;
filtdff = filter(b, a, filt2);
plot(filtdff);
% % ylim([-0.1 0.2])

figure(3);
clf
cutOffFreq = 40;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt3 = x;
filtdff = filter(b, a, filt3);
plot(filtdff);
% ylim([-0.1 0.2])

figure(4);
clf
cutOffFreq = 50;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt4 = x;
filtdff = filter(b, a, filt4);
plot(filtdff);
% ylim([-0.1 0.2])

figure(5);
clf
cutOffFreq = 60;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt5 = x;
filtdff = filter(b, a, filt5);
plot(filtdff);
% ylim([-0.1 0.2])

figure(6);
clf
cutOffFreq = 70;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt6 = x;
filtdff = filter(b, a, filt6);
plot(filtdff);
% ylim([-0.1 0.2])

figure(7)
clf
cutOffFreq = 80;   % Cutoff frequency
[b, a] = butter(4, cutOffFreq/(imagingFreq/2)); % 4th order Butterworth filter
filt7 = x;
filtdff = filter(b, a, filt7);
plot(filtdff);
% ylim([-0.1 0.2])

figure(8);
clf
plot(x);
% ylim([-0.1 0.2])