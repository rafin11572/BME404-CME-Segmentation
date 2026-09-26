function out = segmentCME(img, cfg, stack)
%SEGMENTCME  Run the complete CME segmentation pipeline on one B-scan.
%
%   out = segmentCME(img)              uses the modified (our) configuration
%   out = segmentCME(img, cfg)         uses the supplied configuration
%
% This function is the single entry point used by EVERY other part of the
% project -- the batch scripts, the figure scripts and the GUI all call this.
% No pipeline logic is duplicated anywhere else.
%
% PIPELINE (matching Fig. 1 of Liu et al. 2021)
%     input -> denoise -> contrast stretch     ("working environment")
%           -> ROI extraction                  (retina only)
%           -> omnidirectional wave operator   (4 directional contours)
%           -> contour integration + screening (fused, filtered regions)
%           -> clinical quantification         (OUR addition)
%
% INPUTS
%   img  MxN numeric B-scan.  Accepts uint8 or double; anything whose maximum
%        is <= 1 is assumed to be on [0,1] and is rescaled to 0..255.
%   cfg  config struct from cmeConfig (defaults to cmeConfig('modified'))
%
% OUTPUT (struct)
%   out.mask        MxN logical, the final CME segmentation
%   out.stages      preprocessing intermediates (.raw .denoised .enhanced)
%   out.work        the working image actually fed to the operator
%   out.roi         struct from extractROI
%   out.dirs        struct from omniWaveContours (per-direction + fused)
%   out.integrate   struct from integrateContours (filters + component table)
%   out.clinical    struct from clinicalQuantify (area, count, diameter, grade)
%   out.timing      per-stage seconds, and .total
%   out.cfg         the configuration that produced this result

if nargin < 2 || isempty(cfg), cfg = cmeConfig('modified'); end
if nargin < 3, stack = []; end

I = double(img);
if max(I(:)) <= 1.0 + eps, I = I * 255; end     % accept [0,1] input too

tAll = tic;

% --- Stage 1: construct the working environment ----------------------
t = tic;
[work, stages] = preprocessOCT(I, cfg, stack);
timing.preprocess = toc(t);

% --- Stage 2: ROI extraction -----------------------------------------
t = tic;
roi = extractROI(work, cfg);
timing.roi = toc(t);

% --- Stage 3: omnidirectional wave operator --------------------------
t = tic;
dirs = omniWaveContours(work, cfg, roi.mask);
timing.wave = toc(t);

% --- Stage 4: contour integration and screening ----------------------
t = tic;
% NOTE which image the screening uses.  The mean-gray test of Sec. 2.4 asks
% "is this component hypo-reflective?".  That question must be asked of an
% image whose gray values still mean something physically -- i.e. AFTER
% denoising but BEFORE contrast enhancement.  CLAHE deliberately destroys
% global gray relationships (that is the whole point of local equalisation),
% so screening on the CLAHE output would compare cysts against a rescaled,
% tile-dependent background.
[mask, integrate] = integrateContours(dirs, stages.denoised, roi, cfg, work);
timing.integrate = toc(t);

% --- Stage 5: clinical output layer (ours) ---------------------------
t = tic;
if cfg.clinical.enable
    clinical = clinicalQuantify(mask, roi, cfg);
else
    clinical = struct();
end
timing.clinical = toc(t);

timing.total = toc(tAll);

out.mask      = mask;
out.stages    = stages;
out.work      = work;
out.roi       = roi;
out.dirs      = dirs;
out.integrate = integrate;
out.clinical  = clinical;
out.timing    = timing;
out.cfg       = cfg;
end
