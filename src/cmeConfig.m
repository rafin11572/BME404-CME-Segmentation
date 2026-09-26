function cfg = cmeConfig(mode, varargin)
%CMECONFIG  Build the configuration struct that drives the whole CME pipeline.
%
%   cfg = cmeConfig('baseline')  -> faithful re-implementation of
%                                   Liu et al. (2021), Appl. Sci. 11(14):6480
%   cfg = cmeConfig('modified')  -> our pipeline (median filter + CLAHE +
%                                   adaptive threshold + clinical layer)
%   cfg = cmeConfig(mode,'a.b',value,...) overrides any (possibly nested) field.
%
% INPUTS
%   mode      'baseline' | 'modified'
%   varargin  Name/Value overrides, e.g. cmeConfig('modified','wave.threshold','niblack')
%
% OUTPUT
%   cfg       struct consumed by segmentCME and every module it calls.
%
% DESIGN NOTE
%   Every stage of the pipeline reads its parameters from this one struct, so
%   "baseline mode" and "modified mode" are the *same code path* with
%   different settings.  That is what makes our modifications independently
%   toggleable and separately explainable.
%
% PROVENANCE OF NUMBERS
%   [PAPER] taken directly from Liu et al. (2021) or Lou et al. (2020)
%   [DUKE]  from the Chiu et al. (2015) dataset documentation
%   [LIT]   from the wider clinical literature (cited inline)
%   [OURS]  a default we chose ourselves -- tunable, not claimed from a paper

if nargin < 1 || isempty(mode), mode = 'modified'; end
mode = validatestring(mode, {'baseline','modified'});
cfg.mode = mode;

% ---------------------------------------------------------------------
% 1. Acquisition geometry (Duke / Spectralis 61-line volume protocol)
% ---------------------------------------------------------------------
% Liu et al. (2021) Sec. 3 and Chiu et al. (2015) both state these figures.
cfg.pixelPitchAxialUm   = 3.87;   % [DUKE] microns per pixel, axial   (rows)
cfg.pixelPitchLateralUm = 11.46;  % [DUKE] microns per pixel, lateral (cols)
                                  %        paper quotes 10.94-11.98; midpoint
cfg.pixelAreaMm2 = (cfg.pixelPitchAxialUm*1e-3) * (cfg.pixelPitchLateralUm*1e-3);

% ---------------------------------------------------------------------
% 2. Stage 1 -- construct the working environment (denoise + contrast)
% ---------------------------------------------------------------------
switch mode
    case 'baseline'
        % [PAPER] Liu 2021 Sec. 2.2: "the algorithm in this paper uses
        % Gaussian blur", then "using gamma transformation can stretch the
        % contrast between the tissue and the background".
        cfg.denoise.method   = 'gaussian';
        cfg.denoise.sigma    = 1.0;   % [OURS] paper does not state sigma
        cfg.contrast.method  = 'gamma';
        cfg.contrast.gamma   = 0.8;   % [OURS] <1 brightens; paper only says
                                      %        the stretch must be "moderate"
    case 'modified'
        % [OURS] Proposal modifications 1 and 2.
        cfg.denoise.method   = 'median';
        cfg.denoise.window   = [5 5];
        cfg.contrast.method  = 'clahe';
        cfg.contrast.claheClipLimit = 0.01;  % [OURS] conservative: CLAHE on
                                             %  NOTE the tile count matters a
                                             %  lot: with 8x8 tiles a tile is
                                             %  the same size as a cyst, so
                                             %  CLAHE equalises the cyst away.
                                             %  4x4 keeps tiles safely larger
                                             %  than any lesion.
        cfg.contrast.claheNumTiles  = [4 4]; %        OCT amplifies speckle if
                                             %        the clip limit is high
end
% Both modes keep the same (unset) fields defined so the struct shape is
% identical in either mode -- simplifies the GUI and the comparison scripts.
% Volumetric pre-filter and non-local-means parameters (see preprocessOCT).
% TESTED AND REJECTED -- kept switchable, and the numbers are in
% docs/FINDINGS.md.  The volumetric median genuinely raises the physical edge
% SNR (0.43 -> 0.61) but LOWERS Dice (0.426 -> 0.341), because the ground
% truth is traced on a single slice: anything that mixes information across
% neighbouring B-scans moves our mask away from that slice's tracing.
cfg.denoise.use3D     = false;   % [OURS, tested: worse]
cfg.denoise.shadowComp    = false;  % [OURS, tested: no gain] vessel-shadow fix
cfg.denoise.shadowMaxGain = 1.6;  % [OURS] cap on the per-column brightening
cfg.denoise.nlmDegree = 0.08;   % [OURS] non-local-means smoothing tolerance
cfg.denoise.nlmSearch = 21;     % [OURS] search window side
cfg.denoise.nlmPatch  = 5;      % [OURS] comparison patch side

if ~isfield(cfg.denoise,'sigma'),  cfg.denoise.sigma  = 1.0;   end
if ~isfield(cfg.denoise,'window'), cfg.denoise.window = [5 5]; end
if ~isfield(cfg.contrast,'gamma'), cfg.contrast.gamma = 0.8;   end
if ~isfield(cfg.contrast,'claheClipLimit'), cfg.contrast.claheClipLimit = 0.01; end
if ~isfield(cfg.contrast,"claheNumTiles"),  cfg.contrast.claheNumTiles  = [4 4]; end

% ---------------------------------------------------------------------
% 3. Stage 2 -- ROI extraction (row-mean-gray method)
% ---------------------------------------------------------------------
cfg.roi.smoothRows      = 9;    % [OURS] moving-average span on the row-mean
                                %        curve, so speckle cannot create
                                %        spurious crossings with the mean line
cfg.roi.padTopRows      = 6;    % [OURS] small dilation of the [A,B] band so
cfg.roi.padBottomRows   = 6;    %        we do not clip the ILM / OBM itself
cfg.roi.refineLayers    = true; % refine ILM/OBM per A-scan with the wave op
cfg.roi.layerSmoothCols = 41;   % [OURS] median span for the refined traces
cfg.roi.outlierPx       = 12;   % [OURS] a column whose boundary disagrees with
                                %        its neighbourhood by more than this is
                                %        treated as a detection failure and
                                %        interpolated over
cfg.roi.ilmRiseFrac     = 0.30; % [OURS] the ILM is the first row that rises
                                %        this fraction of the way from the
                                %        VITREOUS background up to the column
                                %        peak.  Referencing the vitreous (not a
                                %        global tissue level) is what stops the
                                %        trace diving into a sub-ILM cyst.
cfg.roi.obmFallFrac     = 0.35; % [OURS] the OBM is where the signal has fallen
                                %        back this fraction below the RPE peak
% TESTED AND REJECTED as a default.  The 'layer' band is the anatomically
% correct one and captures 100% of the expert fluid against 91.2% for
% 'relative' -- but it is 33% larger, and the extra territory is the outer
% nuclear layer, the single worst false-positive region.  Measured over all 78
% B-scans the recall gained is swamped by the precision lost: Dice 0.426
% ('relative') against 0.258 ('layer').  A good lesson that anatomically
% correct is not the same as empirically better.
cfg.roi.bandMode        = 'relative';   % 'layer' | 'relative'  (see extractROI)
cfg.roi.ezOffsetPx      = 11;   % [OURS, measured] the ISM/ISE boundary sits
                                %        this many pixels above the RPE peak
cfg.roi.ilmMarginPx     = 4;    % [OURS] skip the ILM itself
cfg.roi.depthBand       = [0.05 0.55];  % [OURS/PAPER] the "position feature"
                                %        of Liu 2021 Sec 2.4, made explicit:
                                %        CME is confined to this relative depth
                                %        between ILM (0) and OBM (1).  Measured
                                %        on the Duke expert masks: p5 = 0.15,
                                %        median = 0.45, p95 = 0.65.

% ---------------------------------------------------------------------
% 4. Stage 3 -- omnidirectional wave operator
% ---------------------------------------------------------------------
cfg.wave.directionsDeg = [0 90 180 270];  % [PAPER] theta = 0, pi/2, pi, 3pi/2
cfg.wave.gaussSigma    = 0.8;  % [OURS] sigma of the 3x3 Gaussian weighting
                               %        function g(i,j) inside phi_g (Eq. 6).
                               %        The paper names g but never states it.
cfg.wave.epsDiv        = 1e-6; % numerical guard for the divisions in K and C

% How the wave map becomes the binary "wave area" map IT.
%   'alpha'   -> [PAPER] Lou 2020 Eq.(8): IT = wave > alpha, alpha = meanB/meanT
%   'otsu'    -> [OURS]  global Otsu on the wave map
%   'niblack' -> [OURS]  Niblack local threshold T = mean + k*std
switch mode
    case 'baseline', cfg.wave.threshold = 'alpha';
    case 'modified', cfg.wave.threshold = 'otsu';
end
cfg.wave.niblackWindow = 25;    % [OURS] local window side for Niblack
cfg.wave.niblackK      = 0.2;   % [OURS] Niblack k.  Classic Niblack uses
                                %        k = -0.2 for dark text on light paper;
                                %        our feature of interest is BRIGHT on
                                %        the wave map, so k > 0 here.

% Literal ('ratio') vs. sign-safe ('slope') reading of the correction equation.
%   'ratio' -> K = Kb/Kc > 1   exactly as PRINTED in Lou 2020 Eq.(10)
%   'slope' -> Kb > Kc         what the surrounding TEXT says it means
% The two differ whenever Kc < 0.  We default to the literal printed form and
% expose the other so we can show the difference during the viva.
cfg.wave.correctionMode = 'ratio';

% How the boundary line is pulled out of the wave / coastline maps.
%   'edge' : the LEADING EDGE of the wave area, vetted by the correction
%            equation -- a thin line, matching Fig. 5 of the 2021 paper
%   'all'  : every coastline pixel touching the wave area -- the naive literal
%            reading, which marks ~10% of the retina and is unusable
% See the long comment in waveMaps.m for why 'edge' is the right reading.
cfg.wave.boundaryMode = 'edge';

% ---------------------------------------------------------------------
% 5. Stage 4 -- contour integration and screening
% ---------------------------------------------------------------------
% How an open arc becomes a filled region -- the one step the paper leaves
% unspecified.  See the long comment in integrateContours.m.
%   'fill'  : close the arcs + imfill, the paper-literal reading  [BASELINE]
%   'basin' : the cyst is the low-potential-energy BASIN the arcs
%             fence off, recovered by hysteresis                  [MODIFIED]
switch mode
    case 'baseline', cfg.integrate.method = 'fill';
    case 'modified', cfg.integrate.method = 'basin';
end

cfg.integrate.closeRadius   = 8;    % [OURS] only used by the 'fill' strategy.
                                    %        The paper never states how it
                                    %        closes its arcs, so this is a free
                                    %        parameter; 8 px is its best value
                                    %        on our data (3 and 5 close nothing,
                                    %        12 over-merges), and we use the
                                    %        best value so the baseline gets a
                                    %        fair hearing.
cfg.integrate.useFence      = true; % 'basin': let the fused contours separate
                                    %          two cysts that touch
cfg.integrate.fenceDilate   = 1;    % [OURS] thickness of that fence
cfg.integrate.basinImage    = 'denoised'; % 'denoised' | 'enhanced'
cfg.integrate.basinFeature  = 'intensity'; % 'intensity' | 'bothat'
cfg.integrate.bothatShape   = 'disk';
cfg.integrate.bothatSize    = 45;   % [OURS] must exceed the largest cyst
cfg.integrate.thresholdRule = 'sigma';    % 'sigma' | 'percentile' | 'otsu'
cfg.integrate.seedK         = 1.30; % [OURS] seed cut, in sigmas below the mean
cfg.integrate.looseK        = 0.70; % [OURS] loose cut, in sigmas below the mean
cfg.integrate.seedFrac      = 0.55; % [OURS] only used by the 'otsu' rule
cfg.integrate.seedPct       = 4;    % [OURS] strict hysteresis threshold, as a
                                    %        percentile of the gray values in
                                    %        the retina's inner band
cfg.integrate.loosePct      = 20;   % [OURS] loose hysteresis threshold
cfg.integrate.seedMinPx     = 20;   % [OURS] a seed smaller than this is speckle
cfg.integrate.minAreaPx     = 60;   % [OURS] drop specks below this area
cfg.integrate.grayFactor    = 1.0;  % [OURS] multiplier on the paper's mean-gray
                                    %        screening threshold (1.0 = the
                                    %        paper's rule exactly)
cfg.integrate.grayFilter    = true; % [PAPER] Liu 2021 Sec. 2.4 hard-exudate rule
cfg.integrate.thicknessFilter = true;
cfg.integrate.thicknessRange  = [0.5 1.5];  % [PAPER] Liu 2021 Sec. 2.4
cfg.integrate.maxWidthFrac  = 0.60; % [OURS] reject components wider than this
                                    %        fraction of the B-scan: those are
                                    %        dark retinal LAYERS, not cysts
cfg.integrate.maxAspect     = Inf;  % [OURS] max width/height of a component.
                                    %  Cysts are rounded (median w/h 1.6) and
                                    %  dark-layer fragments are long strips
                                    %  (median 2.4), so capping this raises
                                    %  specificity -- but measured over the
                                    %  whole set, maxAspect = 4 halves the
                                    %  false-positive area on healthy scans
                                    %  while costing 17% of Dice (0.445 ->
                                    %  0.370).  We leave it off by default and
                                    %  expose it as an explicit precision /
                                    %  recall trade-off.
cfg.integrate.refineActiveContour = false;  % [OURS] Chan-Vese edge polish
cfg.integrate.acIterations  = 40;   % [OURS] active-contour iterations
cfg.integrate.acSmooth      = 1.0;  % [OURS] contour-length penalty
cfg.integrate.acBias        = 0.0;  % [OURS] >0 shrinks, <0 grows
cfg.integrate.finalErodePx  = 0;    % [OURS] shave the final mask by N px
cfg.integrate.minContrast   = 0;    % [OURS] minimum region-to-surround gray
                                    %        contrast (0 disables the test)

% ---------------------------------------------------------------------
% 6. Stage 5 -- clinical output layer  [OURS, new -- in neither paper]
% ---------------------------------------------------------------------
% See gradeSeverity.m for the full justification of these cut-offs.
cfg.clinical.enable     = true;
cfg.clinical.gradeMode  = 'area';       % 'area' | 'cst' | 'both'
cfg.clinical.areaCutMm2 = [0.05 0.20];  % [OURS, data-driven tertiles]
cfg.clinical.cstCutUm   = [250 400];    % [LIT] DRCR.net / ETDRS-style CST bands
cfg.clinical.centralMm  = 1.0;          % central subfield width used for CST

% ---------------------------------------------------------------------
% 7. Evaluation
% ---------------------------------------------------------------------
if strcmp(mode,'modified')
    % Measured on the Duke expert masks: the median cyst is only ~0.25x the
    % retinal thickness, so the paper's 0.5x lower bound deletes most real
    % cysts on this dataset.  We keep the paper's value in baseline mode and
    % relax it here, and we report the effect as a sensitivity analysis.
    cfg.integrate.thicknessRange = [0.08 1.20];
    cfg.integrate.minContrast    = 10;
    % Chan-Vese polish tested at several contraction biases: neutral to
    % slightly negative, and it costs ~0.5 s per B-scan.  Left off.
    cfg.integrate.refineActiveContour = false;
end

cfg.eval.restrictToROI = true;  % compute TN inside the retina only; scoring
                                % over the whole 496x768 frame inflates
                                % accuracy because most pixels are vitreous.
cfg.eval.restrictToGraded = true;
% The Duke graders annotated only a central window of A-scans (the extent of
% their layer traces).  We verified across all 78 fluid-positive B-scans that
% every annotated fluid pixel lies inside that window, so scoring outside it
% would charge us false positives against ground truth that was never drawn.
% Set false to score over the whole frame instead.

% ---------------------------------------------------------------------
% 8. Apply user overrides -- supports dotted names like 'wave.threshold'
% ---------------------------------------------------------------------
for k = 1:2:numel(varargin)
    parts = strsplit(varargin{k}, '.');
    % subsasgn with a dynamically built substruct lets us poke a nested field
    % without eval().  S is the "subscript reference" describing cfg.a.b.c
    S = struct('type', repmat({'.'},1,numel(parts)), 'subs', parts);
    cfg = subsasgn(cfg, S, varargin{k+1});
end
end
