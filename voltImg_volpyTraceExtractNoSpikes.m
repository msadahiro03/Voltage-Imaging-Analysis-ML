function out = voltImg_volpyTraceExtractNoSpikes(Y, roiMask, varargin)
%VOLTIMG_VOLPYTRACEEXTRACTNOSPIKES VolPy-style preprocessing without spike pursuit.
%   out = voltImg_volpyTraceExtractNoSpikes(Y, roiMask, ...)
%   Y: H x W x T motion-corrected movie (single/double).
%   roiMask: H x W logical, neuron footprint in full-FOV coordinates.
%
%   Pipeline (Cai et al. 2021, PLOS Comp Biol): context crop, optional polarity
%   flip, high-pass (~1/3 Hz) for bleaching, mean ROI trace, SVD of local
%   background + ridge regression to subtract structured background (ridge/SVD in
%   double), one ridge spatial filter over the context (normal equations in double),
%   optional light temporal smoothing, subthreshold
%   via low-pass (default 20 Hz) on the spatially filtered trace.
%
%   Name-value parameters match VolPy defaults where applicable.

p = inputParser;
addParameter(p, 'fps', 400, @isnumeric);
addParameter(p, 'contextDilatePx', 15, @isnumeric);
addParameter(p, 'bgMinDistPx', 12, @isnumeric);  % n_b in VolPy
addParameter(p, 'nPc', 8, @isnumeric);
addParameter(p, 'lambdaB', 0.01, @isnumeric);
addParameter(p, 'lambdaW', 0.01, @isnumeric);
addParameter(p, 'fcHpBleach', 1/3, @isnumeric);
addParameter(p, 'fcLpSubthresh', 20, @isnumeric);
addParameter(p, 'reversePolarity', false, @islogical);
addParameter(p, 'temporalSmoothFrames', 0, @isnumeric); % 0 = off; else movmean width
parse(p, varargin{:});

fps = p.Results.fps;
contextDilatePx = p.Results.contextDilatePx;
bgMinDistPx = p.Results.bgMinDistPx;
nPc = p.Results.nPc;
lambdaB = p.Results.lambdaB;
lambdaW = p.Results.lambdaW;
fcHpBleach = p.Results.fcHpBleach;
fcLpSubthresh = p.Results.fcLpSubthresh;
reversePolarity = p.Results.reversePolarity;
temporalSmoothFrames = p.Results.temporalSmoothFrames;

Y = single(Y);
[H, W, T] = size(Y);
if ~isequal(size(roiMask), [H, W])
    error('voltImg_volpyTraceExtractNoSpikes: roiMask must be H x W.');
end

roiMask = logical(roiMask);
if ~any(roiMask(:))
    error('voltImg_volpyTraceExtractNoSpikes: empty ROI.');
end

ctxMask = imdilate(roiMask, strel('disk', contextDilatePx));
ctxMask = ctxMask & ~isnan(Y(:, :, 1)); % stay in image

bgExclude = imdilate(roiMask, strel('disk', bgMinDistPx));
bgMask = ctxMask & ~bgExclude;
if nnz(bgMask) < max(nPc, 3)
    bgMask = ctxMask & ~roiMask;
end

[rs, cs] = find(ctxMask);
rmin = max(1, min(rs) - 2); rmax = min(H, max(rs) + 2);
cmin = max(1, min(cs) - 2); cmax = min(W, max(cs) + 2);

Yc = Y(rmin:rmax, cmin:cmax, :);
roiC = roiMask(rmin:rmax, cmin:cmax);
ctxC = ctxMask(rmin:rmax, cmin:cmax);
bgC = bgMask(rmin:rmax, cmin:cmax);

if reversePolarity
    Yc = -Yc;
end

Yh = voltImg_temporalFilterVolume(Yc, fps, fcHpBleach, 'high', 3);
[Hc, Wc, ~] = size(Yh);
Nctx = nnz(ctxC);
Nb = nnz(bgC);
if Nctx < 3
    error('voltImg_volpyTraceExtractNoSpikes: context too small.');
end

Ymat = reshape(Yh, Hc * Wc, T)'; % T x (Hc*Wc)
idxCtx = find(ctxC);
idxRoi = find(roiC);
idxBg = find(bgC);

t0 = mean(Ymat(:, idxRoi), 2);
Yb = Ymat(:, idxBg); % T x Nb

nUse = min(nPc, min(size(Yb, 2), T) - 1);
nUse = max(nUse, 1);
if size(Yb, 2) < 2 || T < 4
    tBgSub = t0;
else
    % Double precision: SVD + ridge on background (avoids unstable single solve)
    Yb_d = double(Yb);
    t0_d = double(t0);
    [U, ~, ~] = svd(Yb_d, 'econ');
    U = U(:, 1:min(nUse, size(U, 2)));
    d = size(U, 2);
    A = U' * U + lambdaB * eye(d);
    beta = A \ (U' * t0_d);
    tBgSub = single(t0_d - U * beta);
end

Yctx = Ymat(:, idxCtx); % T x Nctx
d2 = size(Yctx, 2);
% Double precision: normal equations + ridge for spatial weights (large Nctx, single fails)
Yctx_d = double(Yctx);
tBgSub_d = double(tBgSub);
Aw = Yctx_d' * Yctx_d + lambdaW * eye(d2);
bw = Yctx_d' * tBgSub_d;
w = single(Aw \ bw);
tSpatial = Yctx * w;

if temporalSmoothFrames > 1
    k = round(temporalSmoothFrames);
    tSpatial = movmean(tSpatial, k);
end

tSub = voltImg_temporalFilterTrace(tSpatial, fps, fcLpSubthresh, 'low', 5);

out.tMeanRoiHp = t0;
out.tBgSub = tBgSub;
out.tSpatial = tSpatial;
out.tSubthresh = tSub;
out.wSpatial = w;
out.idxCtx = idxCtx;
out.idxRoi = idxRoi;
out.idxBg = idxBg;
out.roiCrop = roiC;
out.ctxCrop = ctxC;
out.bgCrop = bgC;
out.cropBox = [rmin, rmax, cmin, cmax];
end

function y = voltImg_temporalFilterTrace(x, fs, fc, ftype, order)
x = double(x(:));
fn = fc / (fs / 2);
fn = min(max(fn, 1e-6), 0.999999);
[b, a] = butter(order, fn, ftype);
y = single(filtfilt(b, a, x));
end
