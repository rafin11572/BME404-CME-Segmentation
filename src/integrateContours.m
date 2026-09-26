function [mask, info] = integrateContours(dirs, grayImg, roi, cfg, workImg)
%INTEGRATECONTOURS  Stage 4: fuse the directional contours and screen them
%                   into the final CME regions.
%
% Liu 2021 Sec. 2.4 specifies:
%   (1) SUPERIMPOSE the four directional results -- "as long as a pixel is
%       identified as a boundary in a certain direction, it will become a
%       point on the final omnidirectional boundary line";
%   (2) screen connected domains by MEAN GRAY, keeping only those below the
%       image mean, which discards hyper-reflective hard exudate;
%   (3) screen by POSITION and size, keeping only components between 0.5x and
%       1.5x the average retinal thickness.
% It does NOT say how an open arc becomes a filled region.  That gap is where
% the two strategies below differ.
%
% STRATEGY 'fill'  -- the paper-literal reading, used by BASELINE mode.
%   Morphologically close the fused arcs, fill the enclosed holes, screen.
%
% STRATEGY 'basin' -- OUR reading, used by MODIFIED mode.
%   The wave equation is wave = phi_v + phi_g.  A cyst is dark (low phi_g) and
%   flat inside (low phi_v), so it is a LOW-potential-energy basin, and its
%   wall is the ridge the operator traces.  We therefore find the basin bodies
%   by hysteresis on the working image inside the retina's inner band, and use
%   the fused contour map as a FENCE that stops neighbouring basins merging.
%   The operator's output is still doing the work -- it supplies the retinal
%   boundaries, the fence, and the enclosure evidence -- but the cyst body
%   comes from the basin rather than from filling an arc.
%
%   WHY WE NEEDED THIS.  We measured the local edge signal-to-noise across
%   true cyst walls on this dataset: the gray step is ~10 levels against a
%   speckle standard deviation of ~15, i.e. SNR ~0.66 (see docs/FINDINGS.md).
%   Below SNR 1 no purely local edge operator can close a contour reliably, so
%   the literal fill strategy cannot work on the Duke DME data however it is
%   tuned.  We report that honestly rather than hiding it.
%
% INPUTS
%   dirs     struct from omniWaveContours (.fused and .dirBnd)
%   grayImg  MxN double, the DENOISED but not contrast-enhanced image -- the
%            gray screening of step (2) is only meaningful on an image whose
%            gray values still have their original physical meaning
%   roi      struct from extractROI (.mask, .innerBand, .thicknessPx)
%   cfg      config struct
%
% OUTPUTS
%   mask   MxN logical, the final CME segmentation
%   info   struct of intermediates for the figures and the GUI's
%          "why was this component rejected" panel

if nargin < 5, workImg = []; end
fused = dirs.fused;
[M,N] = size(grayImg);
info.fence = fused & roi.mask;

switch lower(cfg.integrate.method)

    % =================================================================
    case 'fill'    % ---- paper-literal ----------------------------------
    % =================================================================
        % imclose = dilation followed by erosion with the same structuring
        % element: it bridges gaps narrower than the element without
        % permanently fattening the arcs.  A disk is used because cyst walls
        % curve in every direction and the disk is the isotropic choice.
        closed = imclose(info.fence, strel('disk', cfg.integrate.closeRadius));
        % imfill(...,'holes') floods inward from the image border and keeps
        % whatever the flood could NOT reach, i.e. every enclosed region.
        filled = imfill(closed, 'holes') & ~closed;
        cand = imopen(filled, strel('disk',2)) & roi.mask;

    % =================================================================
    case 'basin'   % ---- ours -------------------------------------------
    % =================================================================
        band = roi.innerBand;
        if nnz(band) < 200
            mask = false(M,N); info = emptyInfo(info, M, N); return
        end

        % Which image the basin threshold reads.  The gray SCREENING below
        % must use the denoised-but-not-enhanced image (its gray values still
        % mean something physically), but the basin THRESHOLD may benefit from
        % the contrast-enhanced one, so it is a separate choice.
        if strcmpi(cfg.integrate.basinImage,'enhanced') && ~isempty(workImg)
            src = workImg;
        else
            src = grayImg;
        end
        % ---- the feature the basin is found in ----------------------
        % 'intensity'  : the cyst is simply DARK.  Global within the band, so
        %                it is confused by the retina's own layered brightness.
        % 'bothat'     : imbothat(I,SE) = imclose(I,SE) - I, the black top-hat.
        %                Closing with a structuring element LARGER than a cyst
        %                fills the cyst in, so the difference is large exactly
        %                where the image is dark RELATIVE TO ITS OWN LOCAL
        %                surroundings.  That removes the layered background and
        %                the lateral brightness drift in one operation, which
        %                is the textbook reason to reach for a top-hat.
        switch lower(cfg.integrate.basinFeature)
            case 'intensity'
                feat = -src;
            case 'bothat'
                feat = imbothat(src, strel(cfg.integrate.bothatShape, cfg.integrate.bothatSize));
            otherwise
                error('integrateContours:feature','Unknown basinFeature');
        end

        g = feat(band);

        % ---- where to put the two hysteresis levels ------------------
        % Note the sign convention: feat is built so that HIGH means "more
        % likely to be fluid" in either feature, so both cuts are upper ones.
        %
        % A fixed PERCENTILE assumes every scan has the same fluid burden,
        % which is badly wrong here -- across the Duke set the fluid share of
        % the inner retina ranges from 0.9% to 17%.  The 'sigma' rule instead
        % places each cut a fixed number of standard deviations away from the
        % band's own mean, so it adapts per scan; 'otsu' lets Otsu's
        % between-class variance criterion choose the loose cut.
        switch lower(cfg.integrate.thresholdRule)
            case 'percentile'
                tSeed  = prctile(g, 100 - cfg.integrate.seedPct);
                tLoose = prctile(g, 100 - cfg.integrate.loosePct);
            case 'sigma'
                mu = mean(g); sd = std(g);
                tSeed  = mu + cfg.integrate.seedK  * sd;
                tLoose = mu + cfg.integrate.looseK * sd;
            case 'otsu'
                lo = min(g); hi = max(g);
                tLoose = lo + graythresh((g-lo)/max(hi-lo,eps)) * (hi-lo);
                tSeed  = tLoose + cfg.integrate.seedFrac * (hi - tLoose);
            otherwise
                error('integrateContours:rule','Unknown thresholdRule');
        end

        % Hysteresis, exactly the idea Canny uses for edge tracking: a STRICT
        % cut gives seeds that are almost certainly fluid, a LOOSE cut gives
        % the plausible extent, and only loose regions containing a seed
        % survive.  One global cut cannot do both jobs at once.
        seeds = bwareaopen((feat > tSeed) & band, cfg.integrate.seedMinPx);
        loose = imfill(imopen((feat > tLoose) & band, strel('disk',2)), 'holes');

        % The fused contour map is used as a fence so that two cysts separated
        % by a thin bright septum are not reconstructed as one blob.
        if cfg.integrate.useFence
            loose = loose & ~imdilate(info.fence, strel('disk', cfg.integrate.fenceDilate));
        end

        % imreconstruct repeatedly dilates the marker and clips it to the
        % mask until nothing changes, so the result is precisely the union of
        % the loose components that contain at least one seed.
        cand = imreconstruct(seeds & loose, loose);
        cand = imfill(cand, 'holes');

    otherwise
        error('integrateContours:method','Unknown method "%s"', cfg.integrate.method);
end

info.candidates = cand;

% ---------------------------------------------------------------------
% Connected-domain screening
% ---------------------------------------------------------------------
[L, nComp] = bwlabel(cand, 8);
stats = regionprops(L, grayImg, 'Area','MeanIntensity','BoundingBox', ...
                                'PixelIdxList','Centroid','Solidity');

% The gray threshold of Sec. 2.4.  The paper says "the average gray value of
% the whole image"; we take it over the retina, because including the black
% vitreous (roughly half the frame) drags the threshold so low that nothing is
% ever rejected.  Flagged as our interpretation.
imgMean = mean(grayImg(roi.mask)) * cfg.integrate.grayFactor;

thickLo = cfg.integrate.thicknessRange(1) * roi.thicknessPx;
thickHi = cfg.integrate.thicknessRange(2) * roi.thicknessPx;

keep = true(nComp,1);
why  = strings(nComp,1);
area = zeros(nComp,1); mg = area; th = area; wd = area; ct = area;

for c = 1:nComp
    area(c) = stats(c).Area;
    mg(c)   = stats(c).MeanIntensity;
    th(c)   = stats(c).BoundingBox(4);   % axial extent  = "thickness"
    wd(c)   = stats(c).BoundingBox(3);   % lateral extent
    why(c)  = "kept";

    R = false(M,N); R(stats(c).PixelIdxList) = true;

    % region-to-surround contrast: a cyst is a dark cavity IN BRIGHTER TISSUE,
    % so the ring around it must be brighter than its interior.  A dark
    % retinal layer fails this test because its surroundings are dark too.
    ring = imdilate(R, strel('disk',7)) & ~imdilate(R, strel('disk',2)) & roi.mask;
    if nnz(ring) >= 20
        ct(c) = mean(grayImg(ring)) - mean(grayImg(R));
    else
        ct(c) = 0;
    end

    if area(c) < cfg.integrate.minAreaPx
        keep(c) = false; why(c) = "too small"; continue
    end
    if cfg.integrate.grayFilter && mg(c) >= imgMean
        % CME is a fluid-filled cavity -> hypo-reflective -> LOW gray.
        % Hard exudate is a hyper-reflective deposit -> HIGH gray.
        keep(c) = false; why(c) = "too bright (exudate / tissue)"; continue
    end
    if cfg.integrate.thicknessFilter && (th(c) < thickLo || th(c) > thickHi)
        keep(c) = false;
        why(c) = sprintf("thickness %.0f px outside [%.0f %.0f]", th(c), thickLo, thickHi);
        continue
    end
    if cfg.integrate.maxWidthFrac < 1 && wd(c) > cfg.integrate.maxWidthFrac * N
        keep(c) = false; why(c) = "spans the frame (a dark layer, not a cyst)"; continue
    end
    if ct(c) < cfg.integrate.minContrast
        keep(c) = false; why(c) = sprintf("surround contrast %.0f too low", ct(c)); continue
    end
    % Elongation.  A cyst is a rounded cavity; a fragment of a dark retinal
    % LAYER is a long thin horizontal strip.  Measured on labelled components,
    % the width/height ratio separates them better than any intensity feature
    % (true cysts median 1.6, false positives median 2.4, p90 5.5).
    if wd(c)/max(th(c),1) > cfg.integrate.maxAspect
        keep(c) = false;
        why(c) = sprintf("too elongated (w/h = %.1f)", wd(c)/max(th(c),1));
        continue
    end
end

mask = false(M,N);
for c = find(keep)'
    mask(stats(c).PixelIdxList) = true;
end

% ---------------------------------------------------------------------
% Optional boundary refinement with an active contour  [OURS]
% ---------------------------------------------------------------------
% The screening above decides WHICH regions are cysts; it does not polish
% WHERE their edges sit, because the hysteresis mask is built from a single
% global cut.  activecontour evolves the mask boundary to minimise the
% Chan-Vese energy, which rewards a partition whose inside and outside are
% each as uniform as possible.  That is exactly the right model here: a cyst
% is a near-uniform dark region inside near-uniform brighter tissue, and
% Chan-Vese needs no gradient at all -- which matters because the gradient is
% what this dataset does not have (edge SNR 0.6).
%
% 'SmoothFactor' penalises contour length, keeping the result from growing
% speckle-shaped fingers; 'ContractionBias' > 0 shrinks, < 0 grows.
if cfg.integrate.refineActiveContour && any(mask(:))
    src2 = grayImg;
    if strcmpi(cfg.integrate.basinImage,'enhanced') && ~isempty(workImg)
        src2 = workImg;
    end
    refined = activecontour(src2/255, mask, cfg.integrate.acIterations, 'Chan-Vese', ...
        'SmoothFactor', cfg.integrate.acSmooth, ...
        'ContractionBias', cfg.integrate.acBias);
    refined = refined & roi.innerBand;
    % Keep only the refined blobs that still overlap a screened region, so the
    % evolution can adjust a boundary but cannot invent a brand-new cyst.
    refined = imreconstruct(mask & refined, refined);
    refined = bwareaopen(refined, cfg.integrate.minAreaPx);
    if any(refined(:))
        info.beforeRefine = mask;
        mask = refined;
    end
end

% Optional final erosion.  Our recall exceeds our precision, i.e. the masks
% run slightly wide of the expert tracing, so shaving the boundary trades a
% little recall for precision.  Exposed as an explicit knob rather than baked
% in, because which way to lean is a clinical choice.
if cfg.integrate.finalErodePx > 0 && any(mask(:))
    mask = imerode(mask, strel('disk', cfg.integrate.finalErodePx));
    mask = bwareaopen(mask, cfg.integrate.minAreaPx);
end

info.imgMean       = imgMean;
info.thickLimitsPx = [thickLo thickHi];
info.afterArea      = rebuild(stats, area >= cfg.integrate.minAreaPx, M, N);
info.afterGray      = rebuild(stats, area >= cfg.integrate.minAreaPx & mg < imgMean, M, N);
info.afterThickness = mask;

if nComp > 0
    info.compStats = table((1:nComp)', area, mg, ct, th, wd, keep, why, ...
        'VariableNames', {'comp','areaPx','meanGray','contrast','thicknessPx','widthPx','kept','reason'});
else
    info.compStats = emptyStats();
end
end


% =========================================================================
function B = rebuild(stats, sel, M, N)
B = false(M,N);
for c = find(sel(:))'
    B(stats(c).PixelIdxList) = true;
end
end

function info = emptyInfo(info, M, N)
info.candidates = false(M,N);
info.afterArea = false(M,N); info.afterGray = false(M,N); info.afterThickness = false(M,N);
info.imgMean = 0; info.thickLimitsPx = [0 0];
info.compStats = emptyStats();
end

function T = emptyStats()
T = table('Size',[0 8], ...
    'VariableTypes',{'double','double','double','double','double','double','logical','string'}, ...
    'VariableNames',{'comp','areaPx','meanGray','contrast','thicknessPx','widthPx','kept','reason'});
end
