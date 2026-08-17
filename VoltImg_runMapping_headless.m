%% VoltImg_runMapping_headless
% Fully-automatic, single-experiment run of the canonical MultiCell mapping
% pipeline. Edit the REQUIRED fields below, then run this script. No prompts:
% rough ROIs come from Cellpose (auto_roi/voltimg_autoRoi_cellpose.m), QC
% overlay + detection report are written to <saveDirectory>/AutoROI/, and the
% post-analysis plotting cells are skipped (run them later, cell-by-cell,
% after loading the saved results).
%
% One-time setup (internet required once, ~2 GB):
%   bash auto_roi/setup_cellpose_env.sh

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, ...
    fullfile(repoRoot, 'auto_roi'), ...
    fullfile(repoRoot, 'archive_MultiCell'));   % trialExcluder lives in archive_MultiCell

autoCfg = struct();
autoCfg.headless = true;

% ---- REQUIRED: experiment inputs -------------------------------------------
autoCfg.ephysFilePath = '';   % folder with the experiment's DAQ ephys .mat (last .mat in folder is loaded)
autoCfg.imgsFilePath  = '';   % folder with the experiment's imaging TIFFs
autoCfg.UpOrDown      = '2';  % '1' upward GEVI, '2' downward GEVI (char!)
autoCfg.ePhysAvail    = 1;    % 1 = ephys readout available, 2 = none

% ---- Output locations (defaults = the script's current hard-coded paths) ----
autoCfg.savePath   = '/Volumes/X10 Pro/MC Imaging Data';  % motion-corrected TIFFs + checkpoints
autoCfg.resultsDir = ['/Users/masatosadahiro/Documents/Data/Voltage Imaging/' ...
    'Voltage Imaging/voltMapping/Analysis Results/Analysis_newMCanddFFcalc'];
autoCfg.normcorrePath = '';   % '' = NoRMCorre assumed already on the MATLAB path

% ---- Automatic rough-ROI detection (Cellpose) --------------------------------
autoCfg.useAutoRoi = true;
autoCfg.roi = struct( ...
    'pythonExe',         fullfile(repoRoot, 'auto_roi', '.venv_cellpose', 'bin', 'python'), ...
    'wrapperPath',       fullfile(repoRoot, 'auto_roi', 'cellpose_wrapper.py'), ...
    'model',             'cpsam', ...  % Cellpose 4 default (diameter-robust)
    'diameter',          30, ...       % approx soma diameter, px
    'flowThreshold',     0.4, ...
    'cellprobThreshold', 0.0, ...
    'useGpu',            false, ...
    'minAreaPx',         150, ...      % reject tight masks smaller than this
    'maxAreaPx',         5000, ...     % reject merged blobs / debris larger than this
    'minSeparationPx',   15, ...       % greedy: drop dimmer cell of any closer pair
    'maxCells',          Inf, ...      % optional cap (keeps brightest N)
    'excludeBorderPx',   0, ...        % reject cells touching an N-px border margin (0 = off)
    'dilateRadiusPx',    3);           % expand tight masks to emulate hand-drawn margins (3 px best matched hand ROIs in validation)

if isempty(autoCfg.ephysFilePath) || isempty(autoCfg.imgsFilePath)
    error('VoltImg:headless', 'Set autoCfg.ephysFilePath and autoCfg.imgsFilePath first.');
end

VoltImg_mapping_analysis_MultiCell_newDFF_021226_MCfineROI_TrialSpecROIandNeuropil_laserRowArtifact
