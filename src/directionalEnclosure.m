function [score, perDir] = directionalEnclosure(regionMask, dirBnd, roiMask, maxGapPx)
%DIRECTIONALENCLOSURE  How completely is a candidate region walled in, in each
%                      of the operator's four directions?
%
% THIS IS THE STEP THAT MAKES THE OMNIDIRECTIONAL OPERATOR PAY FOR ITSELF.
%
% Liu et al. (2021) argue that four directions are needed because a closed
% curve cannot be traced from one direction alone.  The converse is the useful
% test: a genuine cyst is a closed curve, so it must have a wave boundary on
% ALL FOUR sides.  A dark retinal LAYER (the outer nuclear layer is genuinely
% hypo-reflective and is the main false positive on this dataset) has walls
% above and below but none to the left or right -- it runs off the edge of the
% frame.  Scoring enclosure per direction separates the two cleanly, using
% nothing but the directional contour maps the operator already produced.
%
% INPUTS
%   regionMask  MxN logical, ONE candidate region
%   dirBnd      1x4 cell of MxN logical directional boundary maps, in the
%               order [0 90 180 270] degrees as produced by omniWaveContours
%   roiMask     MxN logical retina mask
%   maxGapPx    how far to look for a wall beyond the region edge (default 12)
%
% OUTPUTS
%   perDir  1x4, the fraction of the region's exit points in each direction
%           that have a wall within maxGapPx
%   score   min(perDir) -- a region is only "closed" if its WEAKEST side is
%           walled, so the minimum is the right summary, not the mean
%
% Direction convention (matching omniWaveContours / waveMaps):
%   theta =   0  -> operator sweeps DOWN  (+row); its contour is the FLOOR
%   theta =  90  -> sweeps one way along the columns; contour is a SIDE wall
%   theta = 180  -> sweeps UP    (-row); its contour is the ROOF
%   theta = 270  -> the other SIDE wall

if nargin < 4 || isempty(maxGapPx), maxGapPx = 12; end

% Each direction is checked by walking from the region's exit pixels outwards.
% Rotating the problem (rather than writing four near-identical loops) keeps
% this short and guarantees the four tests are identical in everything but
% direction -- exactly the same trick omniWaveContours uses.
perDir = zeros(1,4);
for t = 1:4
    k = t - 1;                       % number of 90-degree turns
    R  = rot90(regionMask, k);
    Bd = rot90(dirBnd{t},  k);
    Ro = rot90(roiMask,    k);
    perDir(t) = sideSupport(R, Bd, Ro, maxGapPx);
end
score = min(perDir);
end


% =========================================================================
function frac = sideSupport(R, B, roi, maxGap)
%SIDESUPPORT  Fraction of the region's LOWER edge that has a wall below it.
% After the rotation above, every direction reduces to this one case.
[M,~] = size(R);
cols = find(any(R,1));
if isempty(cols), frac = 0; return; end

hit = 0; tot = 0;
for c = cols
    r = find(R(:,c), 1, 'last');     % the region's exit point in this column
    if isempty(r), continue; end
    tot = tot + 1;
    lo = r + 1; hi = min(M, r + maxGap);
    if lo > M
        % The region runs into the image border: treat the border as a wall,
        % otherwise every region touching the frame edge is scored as open.
        hit = hit + 1; continue
    end
    seg = B(lo:hi, c);
    if any(seg)
        hit = hit + 1;
    elseif ~any(roi(lo:hi, c))
        % We left the retina before finding a wall -- the ROI edge is itself a
        % legitimate boundary (a cyst can sit right under the ILM).
        hit = hit + 1;
    end
end
frac = hit / max(tot,1);
end
