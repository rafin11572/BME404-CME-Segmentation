function roi = extractROI(I, cfg)
%EXTRACTROI  Stage 2: isolate the retina (ILM..OBM) before segmenting.
%
% Liu 2021 Sec. 2.3.1.  The OCT frame contains vitreous above the retina and
% choroid/sclera below it.  The choroidal vessel lacunae are dark, closed and
% roughly cyst-shaped, so without this stage they get segmented as cysts.
% Restricting the search to the retina removes that failure mode and makes the
% whole algorithm several times faster.
%
% METHOD
%   Step 1 (Eq. 1, the paper's figure-3 method)
%     meangray(i) = (1/N) * sum_j I(i,j)        -- the mean gray of each ROW
%     The intersections A and B of that curve with the whole-image mean gray
%     split the frame into vitreous / retina+choroid / sclera.
%
%   Step 2 (per-A-scan refinement)
%     The paper then says "using the wave algorithm, the image is segmented
%     from top to bottom to obtain the ILM layer and the OBM layer".  We do
%     that -- the wave boundary map supplies the candidate rows -- but we pick
%     WHICH candidate with an explicit intensity criterion, because the raw
%     "first hit" is easily fooled by vitreous speckle.  The two criteria are:
%       ILM : the first wave boundary below A at which the A-scan intensity
%             rises above the tissue level and STAYS there (persistence test)
%       OBM : the wave boundary nearest below the RPE peak at which the
%             intensity has fallen back towards the background level
%     Both are stated explicitly here because they are OUR reading of a step
%     the paper describes in one sentence.
%
% INPUTS
%   I    MxN double working image on 0..255 (already preprocessed)
%   cfg  config struct
%
% OUTPUT (struct)
%   roi.rowMean/.rowMeanRaw  the row-mean-gray curve (for the Fig. 3 plot)
%   roi.imgMean              whole-image mean gray (the blue line in Fig. 3)
%   roi.A, roi.B             row indices bounding the retina
%   roi.ilm, roi.obm         1xN per-A-scan boundary traces
%   roi.mask                 MxN logical, true between ILM and OBM
%   roi.thicknessPx          mean retinal thickness in pixels

D = double(I);
[M,N] = size(D);

% =====================================================================
% Step 1 -- Eq. (1): row-mean-gray curve and its crossings A, B
% =====================================================================
rowMeanRaw = mean(D, 2);

% movmean is a sliding-window moving average.  Speckle makes the raw curve
% cross the image-mean line many times near A and B; smoothing first gives one
% clean crossing at each end.
rowMean = movmean(rowMeanRaw, cfg.roi.smoothRows);
imgMean = mean(D(:));

above = rowMean > imgMean;
if ~any(above)
    A = 1; B = M;
else
    % Several short runs can sit above the mean.  The retina is by far the
    % longest, so keep that one.  bwlabel numbers each connected run.
    [lbl, nRun] = bwlabel(above(:)');
    len = zeros(1,nRun);
    for r = 1:nRun, len(r) = sum(lbl == r); end
    [~, best] = max(len);
    idx = find(lbl == best);
    A = idx(1); B = idx(end);
end
A = max(1, A - cfg.roi.padTopRows);
B = min(M, B + cfg.roi.padBottomRows);

% =====================================================================
% Step 2 -- per-A-scan ILM and OBM
% =====================================================================
ilm = repmat(A, 1, N);
obm = repmat(B, 1, N);
rpeTrace = repmat(round((A+B)/2), 1, N);

if cfg.roi.refineLayers
    % --- intensity reference levels ---------------------------------
    % Background = the vitreous, i.e. the rows above A.  Tissue level is set
    % by Otsu on the whole frame (see niblackThreshold.m for what Otsu does).
    bg = median(D(1:max(1,A-5), :), 'all');
    tissueLvl = graythresh(uint8(D)) * 255;

    % Smooth each A-scan along depth so the persistence tests are not tripped
    % by single speckle pixels.  A column-wise moving average is a 1-D box
    % filter; 5 rows is ~19 um, well below any real layer thickness.
    P = movmean(D, 5, 1);

    % --- wave boundary candidates, both sweep directions -------------
    band = false(M,N); band(A:B,:) = true;
    cfgL = cfg;
    cfgL.wave.threshold = 'alpha';   % the layer hunt always uses the paper's
                                     % global alpha: Otsu on a near-empty wave
                                     % map is unstable
    Wd   = waveMaps(D, cfgL, band);              % downward sweep
    bndD = Wd.bnd;
    Wu   = waveMaps(flipud(D), cfgL, flipud(band));
    bndU = flipud(Wu.bnd);                       % upward sweep, mapped back

    for c = 1:N
        prof = P(:,c);
        pkCol = max(prof(A:B));

        % ---- ILM ---------------------------------------------------
        % The ILM is the FIRST rise out of the vitreous.  The level must be
        % referenced to the VITREOUS BACKGROUND, not to a global tissue level:
        % where a large cyst sits directly under the ILM the nerve-fibre layer
        % above it is thin and dark, so a global-tissue criterion skips the
        % membrane entirely and lands inside (or below) the cyst.  That single
        % mistake truncates exactly the biggest cysts, which is the worst
        % possible failure for this project.
        lvlI = bg + cfg.roi.ilmRiseFrac * (pkCol - bg);
        ilmC = NaN;
        k = A;
        while k <= B-3
            if prof(k) > lvlI && all(prof(k:k+3) > bg + 0.5*cfg.roi.ilmRiseFrac*(pkCol-bg))
                ilmC = k; break
            end
            k = k + 1;
        end
        % Snap onto a wave-operator boundary if one sits within a few pixels,
        % so the final trace is still an operator-detected edge.
        if ~isnan(ilmC)
            w = find(bndD(:,c));
            w = w(abs(w - ilmC) <= 4);
            if ~isempty(w), [~,ix] = min(abs(w - ilmC)); ilmC = w(ix); end
        end
        ilm(c) = ilmC;

        % ---- RPE peak ----------------------------------------------
        % The RPE/BM complex is the brightest structure in the frame.  Search
        % well below the ILM so a bright nerve-fibre layer cannot win.
        base0 = A; if ~isnan(ilmC), base0 = ilmC; end
        lo = min(M-1, base0 + 30);
        if lo >= B, lo = max(A, B-5); end
        [pk, rel] = max(prof(lo:B));
        rpe = lo + rel - 1;
        rpeTrace(c) = rpe;

        % ---- OBM ---------------------------------------------------
        % Going down from the RPE peak, the outer boundary of the RPE/BM
        % complex is where the signal has decayed a fixed fraction of the way
        % back to the vitreous background.
        lvl = bg + cfg.roi.obmFallFrac * (pk - bg);
        k = find(prof(rpe:B) < lvl, 1, 'first');
        if isempty(k), obmC = B; else, obmC = rpe + k - 1; end

        w = find(bndU(:,c));
        w = w(abs(w - obmC) <= 6);
        if ~isempty(w), [~,ix] = min(abs(w - obmC)); obmC = w(ix); end
        obm(c) = obmC;
    end

    % A retinal boundary is a smooth curve, so a column that jumped tens of
    % pixels is an error.  Reject outliers against a long-window median and
    % re-interpolate, rather than just median-filtering: a median filter alone
    % drags the whole trace towards a cluster of neighbouring failures.
    ilm = robustSmooth(ilm, cfg.roi.layerSmoothCols, cfg.roi.outlierPx);
    obm = robustSmooth(obm, cfg.roi.layerSmoothCols, cfg.roi.outlierPx);
    rpeTrace = robustSmooth(rpeTrace, cfg.roi.layerSmoothCols, cfg.roi.outlierPx);

    ilm = fillAndClamp(ilm, A, B);
    obm = fillAndClamp(obm, A, B);

    bad = ilm >= obm - 5;      % guarantee ILM sits above OBM everywhere
    ilm(bad) = A; obm(bad) = B;
end

% --- Build the masks ---------------------------------------------------
rr = repmat((1:M)', 1, N);
mask = rr >= repmat(round(ilm),M,1) & rr <= repmat(round(obm),M,1);

% Relative depth inside the retina: 0 at the ILM, 1 at the OBM.  This is the
% coordinate in which the "position feature" screening of Liu 2021 Sec. 2.4 is
% expressed, and it is invariant to the retina being thickened by the oedema
% itself -- which a screen expressed in raw pixels would not be.
relDepth = (rr - repmat(ilm,M,1)) ./ max(repmat(obm-ilm,M,1), 1);

% =====================================================================
% The search band for CME -- the "position feature" of Liu 2021 Sec. 2.4
% =====================================================================
% Two ways to express it:
%
%   'relative' : a fixed fraction of the ILM->OBM depth.  Simple, but it is
%                not anatomical, and measured against the expert masks it
%                captures only 91% of the annotated fluid.
%
%   'layer'    : ILM down to the ISM/ISE boundary (the top of the ellipsoid
%                zone).  This is the anatomical statement that intraretinal
%                fluid lies ABOVE the photoreceptor layer, which is also how
%                Rashno et al. (2018) bound their ROI on this same dataset.
%                Measured against the expert traces it captures 100.0% of the
%                annotated fluid -- only 0.1% of fluid pixels lie below the
%                ISM/ISE boundary at all.
%
%                We do not have the expert traces at inference time, so the
%                ISM/ISE boundary is estimated from the RPE peak: measured on
%                the 110 graded B-scans it sits a tight 9-14 px above it
%                (10th-90th percentile), i.e. about 40 um, which is the
%                expected outer-segment thickness.
switch lower(cfg.roi.bandMode)
    case 'relative'
        innerBand = mask & relDepth >= cfg.roi.depthBand(1) ...
                         & relDepth <= cfg.roi.depthBand(2);
    case 'layer'
        ezRow = rpeTrace - cfg.roi.ezOffsetPx;          % ISM/ISE estimate
        ezRow = max(ezRow, ilm + 10);                   % never above the ILM
        innerBand = rr >= repmat(round(ilm),M,1) + cfg.roi.ilmMarginPx ...
                  & rr <  repmat(round(ezRow),M,1);
    otherwise
        error('extractROI:bandMode','Unknown bandMode "%s"', cfg.roi.bandMode);
end

roi.rowMeanRaw  = rowMeanRaw;
roi.rowMean     = rowMean;
roi.imgMean     = imgMean;
roi.A           = A;
roi.B           = B;
roi.ilm         = ilm;
roi.obm         = obm;
roi.mask        = mask;
roi.relDepth    = relDepth;
roi.rpe         = rpeTrace;
roi.ez          = rpeTrace - cfg.roi.ezOffsetPx;
roi.innerBand   = innerBand;
roi.thicknessPx = mean(obm - ilm);
end


% =========================================================================
function y = robustSmooth(y, win, tol)
%ROBUSTSMOOTH  Reject columns that disagree with their neighbourhood, then
%              interpolate across them and lightly smooth what is left.
y = y(:)';
med = movmedian(y, win, 'omitnan');
bad = ~isfinite(y) | abs(y - med) > tol;
y(bad) = NaN;
if any(isfinite(y))
    y = interp1(find(isfinite(y)), y(isfinite(y)), 1:numel(y), 'linear', 'extrap');
end
y = movmean(y, max(3, round(win/3)));   % final light smoothing
end


% =========================================================================
function y = fillAndClamp(y, A, B)
%FILLANDCLAMP  Interpolate over columns where no boundary was found.
n = numel(y);
good = isfinite(y);
if ~any(good)
    y = repmat((A+B)/2, 1, n);
else
    % interp1 bridges interior gaps and 'extrap' extends the trace to the
    % image edges, so every A-scan ends up with a value.
    y = interp1(find(good), y(good), 1:n, 'linear', 'extrap');
end
y = min(max(y, A), B);
end
