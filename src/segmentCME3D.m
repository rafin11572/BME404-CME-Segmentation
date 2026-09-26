function out = segmentCME3D(subjectId, k, cfg, opts)
%SEGMENTCME3D  Segment one B-scan using its NEIGHBOURING slices as evidence.
%
%   out = segmentCME3D(5, 24)
%   out = segmentCME3D(5, 24, cmeConfig('modified'), opts)
%
% ---------------------------------------------------------------------
% WHY THIS IS NOT THE SAME AS THE 3-D DENOISING WE ALREADY REJECTED
% ---------------------------------------------------------------------
%   docs/FINDINGS.md records that taking a per-pixel MEDIAN across B-scans
%   k-1, k, k+1 BEFORE segmenting made things worse (Dice 0.426 -> 0.341).
%   That mixes the IMAGES, so a cyst that shifts slightly between slices gets
%   smeared and the mask drifts away from the tracing drawn on slice k.
%
%   This function does something different: it segments each slice
%   INDEPENDENTLY and then uses the neighbouring RESULTS as supporting
%   evidence for the centre slice.  The images are never mixed.  The idea is
%   the one a human grader uses -- a real cyst is a three-dimensional pocket
%   of fluid, so it appears in the slices either side of the one you are
%   looking at, whereas a speckle-driven false positive does not.
%
% ---------------------------------------------------------------------
% INPUTS
%   subjectId  1..10
%   k          B-scan index 1..61
%   cfg        config struct (default cmeConfig('modified'))
%   opts       struct, all optional:
%     .mode          'none'      pass the 2-D result straight through
%                    'support'   keep a component only if enough of it is
%                                backed by a neighbour               (default)
%                    'intersect' keep a PIXEL only if a neighbour has it too
%                    'majority'  keep a pixel present in >= 2 of the 3 slices
%                    'union'     add fluid both neighbours agree on, even if
%                                the centre slice missed it (recall-oriented)
%     .halfWidth     how many slices either side to use          (default 1)
%     .tolerancePx   dilation applied to a neighbour mask before comparing,
%                    to allow for the retina shifting slightly between
%                    slices                                       (default 3)
%     .support       for 'support': the fraction of a component that must be
%                    backed by at least one neighbour             (default 0.30)
%
% OUTPUT
%   out          the same struct segmentCME returns for the centre slice,
%                with .mask replaced by the volumetrically screened mask, plus
%   out.mask2D   the original 2-D mask, so the two can be compared
%   out.support  per-component neighbour support fraction
%   out.neighbourMasks, out.neighbourIdx
%
% COST
%   One segmentCME call per slice used, so roughly 3x the 2-D runtime at the
%   default half-width of 1.  Nothing else changes.

if nargin < 3 || isempty(cfg),  cfg  = cmeConfig('modified'); end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mode'),        opts.mode        = 'support'; end
if ~isfield(opts,'halfWidth'),   opts.halfWidth   = 1;         end
if ~isfield(opts,'tolerancePx'), opts.tolerancePx = 3;         end
if ~isfield(opts,'support'),     opts.support     = 0.30;      end

% ---------------------------------------------------------------------
% 1. Segment the centre slice and each neighbour, independently
% ---------------------------------------------------------------------
Sc  = loadBScan(subjectId, k, 1);
out = segmentCME(Sc.img, cfg, Sc.stack);
out.mask2D = out.mask;
out.scan   = Sc;

nSlices = size(Sc.stack, 3);            % how many exist in the volume
kList   = neighbourIndices(subjectId, k, opts.halfWidth);
kList   = setdiff(kList, k);

nb = {};
for j = 1:numel(kList)
    Sn = loadBScan(subjectId, kList(j), 1);
    rn = segmentCME(Sn.img, cfg, Sn.stack);
    % Restrict to the centre slice's retina: the ROI moves a little between
    % slices and we only ever ask "is there fluid HERE", in centre-slice
    % coordinates.
    nb{end+1} = rn.mask; %#ok<AGROW>
end
out.neighbourMasks = nb;
out.neighbourIdx   = kList;

if isempty(nb) || strcmpi(opts.mode,'none')
    out.support = [];
    return
end

% ---------------------------------------------------------------------
% 2. Build the neighbour evidence map
% ---------------------------------------------------------------------
% Dilate each neighbour mask slightly before comparing.  Consecutive B-scans
% are ~118 um apart azimuthally, so the same cyst is not at pixel-identical
% coordinates in the next slice; without a tolerance the comparison would be
% far too strict and would delete genuine fluid.
se = strel('disk', opts.tolerancePx);
votes = zeros(size(out.mask));                       % how many slices agree
for j = 1:numel(nb)
    votes = votes + double(imdilate(nb{j}, se));
end
out.votes = votes;

% ---------------------------------------------------------------------
% 3. Apply the chosen rule
% ---------------------------------------------------------------------
m2 = out.mask2D;
switch lower(opts.mode)

    case 'support'
        % Component-level.  A component survives only if a decent fraction of
        % it is backed by at least one neighbour.  Working per component
        % rather than per pixel means a real cyst keeps its full extent --
        % we are deciding "is this thing real", not "which of its pixels are".
        [L, nComp] = bwlabel(m2, 8);
        st  = regionprops(L, 'PixelIdxList', 'Area');
        keep = false(size(m2));
        sup  = zeros(nComp,1);
        for c = 1:nComp
            idx    = st(c).PixelIdxList;
            sup(c) = nnz(votes(idx) > 0) / numel(idx);
            if sup(c) >= opts.support
                keep(idx) = true;
            end
        end
        out.mask    = keep;
        out.support = sup;

    case 'intersect'
        % Pixel-level, strict: the centre slice AND at least one neighbour.
        out.mask = m2 & (votes > 0);

    case 'majority'
        % Pixel-level: present in at least two of the three slices.  The
        % centre slice counts as one vote.
        out.mask = (double(m2) + votes) >= 2;

    case 'union'
        % Recall-oriented: keep everything the centre slice found, and ADD
        % anything that BOTH neighbours agree on, on the grounds that fluid
        % present either side is very likely present in between.  Restricted
        % to the centre slice's own search band so it cannot invent fluid
        % outside the retina.
        addIn    = (votes >= numel(nb)) & out.roi.innerBand;
        out.mask = m2 | addIn;

    otherwise
        error('segmentCME3D:mode','Unknown mode "%s"', opts.mode);
end

% Re-apply the size floor so no rule can leave slivers behind.
out.mask = bwareaopen(out.mask, cfg.integrate.minAreaPx);

% The clinical numbers must describe the mask we actually return.
if cfg.clinical.enable
    out.clinical = clinicalQuantify(out.mask, out.roi, cfg);
end
end


% =========================================================================
function kk = neighbourIndices(subjectId, k, halfWidth)
%NEIGHBOURINDICES  Slice indices around k that exist in this volume.
p = dukePaths();
f = fullfile(p.data, 'v73', sprintf('Subject_%02d.mat', subjectId));
if exist(f,'file') ~= 2
    f = fullfile(p.data, sprintf('Subject_%02d.mat', subjectId));
end
m  = matfile(f);
nK = size(m, 'images', 3);
kk = max(1, k-halfWidth) : min(nK, k+halfWidth);
end
