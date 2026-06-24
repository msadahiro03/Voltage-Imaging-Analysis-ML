function Yf = voltImg_temporalFilterVolume(Y, fs, fc, ftype, order)
%VOLTIMG_TEMPORALFILTERVOLUME Apply the same temporal IIR filter to every pixel.
%   Y: H x W x T, fs in Hz, fc cutoff (Hz), ftype 'high'|'low', order default 3.

if nargin < 5
    order = 3;
end
Y = single(Y);
[H, W, T] = size(Y);
fn = fc / (fs / 2);
fn = min(max(fn, 1e-6), 0.999999);
[b, a] = butter(order, fn, ftype);
np = H * W;
Y2 = reshape(permute(Y, [3 1 2]), T, np); % T x np
Y2f = zeros(T, np, 'single');
for k = 1:np
    Y2f(:, k) = single(filtfilt(b, a, double(Y2(:, k))));
end
Yf = permute(reshape(Y2f, T, H, W), [2 3 1]);
end
