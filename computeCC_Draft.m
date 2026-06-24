function out = compute_RDI(session1_trials_by_frames, sessionN_trials_by_frames, n_iter)
% sessionX_trials_by_frames: [nTrials x nFrames] ΔF/F (or inferred spikes) for ONE neuron
% n_iter: number of random half-split resamples (e.g., 100–500)
%
% Returns:
%   out.ccws_s1   - within-session CC for Session 1 (averaged over n_iter)
%   out.ccws_sn   - within-session CC for Session N (averaged over n_iter)
%   out.ccbs      - between-session CC (averaged over n_iter)
%   out.rdi       - Representational Drift Index
%
% Method matches the paper:
% - CCws: Pearson correlation between trial-averaged activity of two random halves WITHIN a session
% - CCbs: Pearson correlation between trial-averaged activity of two random halves ACROSS the two sessions
% - Negative correlations are rectified to 0 before use
% - RDI = (CCws - CCbs) / (CCws + CCbs)  [using CCws from the reference session]
%
% Refs:
%   RDI / CCws / CCbs definitions: Marks & Goard 2021, Methods. 
%   (trial-averaged activity, halves, rectification, and RDI formula)

if nargin < 3, n_iter = 200; end

S1 = session1_trials_by_frames;
SN = sessionN_trials_by_frames;

assert(size(S1,2) == size(SN,2), 'Sessions must have same number of frames.');

% convenience
n1 = size(S1,1);
nN = size(SN,1);
half1 = floor(n1/2);
halfN = floor(nN/2);

ccws_s1_vals = nan(n_iter,1);
ccws_sn_vals = nan(n_iter,1);
ccbs_vals    = nan(n_iter,1);

for k = 1:n_iter
    % ---- Session 1: CCws (within) ----
    idx = randperm(n1);
    a = idx(1:half1);
    b = idx(half1+1 : 2*half1);  % if odd trials, 1 trial is unused
    mA = mean(S1(a,:), 1);
    mB = mean(S1(b,:), 1);
    c = corr_no_nan(mA(:), mB(:));
    if c < 0, c = 0; end
    ccws_s1_vals(k) = c;

    % ---- Session N: CCws (within) ----
    idx = randperm(nN);
    a = idx(1:halfN);
    b = idx(halfN+1 : 2*halfN);
    mA = mean(SN(a,:), 1);
    mB = mean(SN(b,:), 1);
    c = corr_no_nan(mA(:), mB(:));
    if c < 0, c = 0; end
    ccws_sn_vals(k) = c;

    % ---- CCbs (between sessions) ----
    % Split each session independently into halves, then correlate across matching halves.
    % Average the two possible pairings to be symmetric.
    % (You can also use one pairing; averaging is a small-variance nicety.)
    % S1 halves:
    idx1 = randperm(n1);
    s1a = idx1(1:half1);  s1b = idx1(half1+1 : 2*half1);
    m1a = mean(S1(s1a,:), 1);  m1b = mean(S1(s1b,:), 1);
    % SN halves:
    idxN = randperm(nN);
    sna = idxN(1:halfN);  snb = idxN(halfN+1 : 2*halfN);
    mNa = mean(SN(sna,:), 1);  mNb = mean(SN(snb,:), 1);

    c1 = corr_no_nan(m1a(:), mNa(:));
    c2 = corr_no_nan(m1b(:), mNb(:));
    c_pair = mean([max(c1,0), max(c2,0)]);  % rec_
