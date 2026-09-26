function R = omniWaveContours(I, cfg, analysisMask)
%OMNIWAVECONTOURS  Run the wave operator in every operating direction and
%                  return the individual and fused contour maps.
%
% This is the "direction adjustment function" of Liu 2021 Eq. (5):
%
%     [i']   [ cos(th)  -sin(th)  0 ] [i]
%     [j'] = [ sin(th)   cos(th)  0 ] [j]
%     [1 ]   [   0         0      1 ] [1]
%
% Rotating the COORDINATES by theta is the same as rotating the IMAGE by
% -theta and then running the fixed, canonical (downward) operator on it.  We
% do the latter because it lets waveMaps.m stay simple and readable: it only
% ever has to know about one direction.
%
% INPUTS
%   I             MxN double working image, 0..255
%   cfg           config struct (cfg.wave.directionsDeg lists the angles)
%   analysisMask  MxN logical ROI; contours outside it are discarded
%
% OUTPUTS (struct)
%   R.dirBnd    1xD cell, each MxN logical: the contour found in one direction
%               (these are Fig. 5(a)-(d) of the 2021 paper)
%   R.fused     MxN logical: the union over directions
%               (Fig. 5(e) -- "as long as a pixel is identified as a boundary
%                in a certain direction, it will become a point on the final
%                omnidirectional boundary line")
%   R.anglesDeg the angles actually used
%   R.exemplar  the full waveMaps struct for the theta = 0 sweep, kept so the
%               figure scripts can draw the seawater/wave/coastline maps
%
% WHY FOUR DIRECTIONS AND NOT SIX OR EIGHT
%   For theta a multiple of 90 degrees the rotation is a pure re-indexing
%   (rot90) -- no pixel is invented and no gray value is interpolated.  For 6
%   or 8 directions the rotation needs interpolation, which, as the paper puts
%   it, "will change the image pixel distribution to a certain extent and
%   affect the segmentation result".  Our implementation makes that concrete:
%   exact rot90 below vs. bilinear imrotate.  This is why the ablation in
%   run_04 is not just a number, it is a statement about image authenticity.

if nargin < 3 || isempty(analysisMask), analysisMask = true(size(I)); end

ang = cfg.wave.directionsDeg;
R.anglesDeg = ang;
R.dirBnd = cell(1, numel(ang));
fused = false(size(I));

for t = 1:numel(ang)
    th = ang(t);

    % --- rotate into the operator's canonical frame -------------------
    [Ir, Mr, meta] = rotateToCanonical(I, analysisMask, th);

    % --- run the fixed downward-sweeping operator ---------------------
    W = waveMaps(Ir, cfg, Mr);

    % --- rotate the boundary map back into image coordinates ----------
    bnd = rotateFromCanonical(W.bnd, meta);
    bnd = bnd & analysisMask;

    R.dirBnd{t} = bnd;
    fused = fused | bnd;          % union across directions

    if th == 0, R.exemplar = W; end
end

if ~isfield(R,'exemplar')
    R.exemplar = waveMaps(I, cfg, analysisMask);
end
R.fused = fused;
end


% =========================================================================
function [Ir, Mr, meta] = rotateToCanonical(I, M, thetaDeg)
%ROTATETOCANONICAL  Put the operating direction for theta onto the +row axis.
%
% For multiples of 90 degrees we use rot90, which just permutes existing
% pixels -- lossless, and the reason the paper prefers four directions.
% For any other angle we pad the image to a square big enough that nothing
% rotates out of frame, then use bilinear imrotate.

meta.thetaDeg = thetaDeg;
meta.origSize = size(I);
q = thetaDeg / 90;

if abs(q - round(q)) < 1e-9
    meta.kind = 'rot90';
    meta.k    = mod(round(q), 4);
    Ir = rot90(I, meta.k);
    Mr = rot90(M, meta.k);
else
    meta.kind = 'imrotate';
    [m,n] = size(I);
    s = ceil(hypot(m,n));                 % square side that can hold any rotation
    meta.pad = [floor((s-m)/2), floor((s-n)/2)];
    meta.sq  = s;
    Ip = zeros(s); Mp = false(s);
    r0 = meta.pad(1)+1; c0 = meta.pad(2)+1;
    Ip(r0:r0+m-1, c0:c0+n-1) = I;
    Mp(r0:r0+m-1, c0:c0+n-1) = M;
    % imrotate rotates counter-clockwise by the given angle. 'bilinear'
    % interpolates each output pixel from the 4 nearest input pixels;
    % 'crop' keeps the output the same size as the (already oversized) input.
    Ir = imrotate(Ip, thetaDeg, 'bilinear', 'crop');
    Mr = imrotate(Mp, thetaDeg, 'nearest',  'crop') > 0;  % nearest keeps the
                                                          % mask strictly binary
end
end


% =========================================================================
function B = rotateFromCanonical(Br, meta)
%ROTATEFROMCANONICAL  Undo rotateToCanonical for a binary map.
switch meta.kind
    case 'rot90'
        B = rot90(Br, -meta.k);
    case 'imrotate'
        Bp = imrotate(double(Br), -meta.thetaDeg, 'nearest', 'crop') > 0;
        m = meta.origSize(1); n = meta.origSize(2);
        r0 = meta.pad(1)+1; c0 = meta.pad(2)+1;
        B = Bp(r0:r0+m-1, c0:c0+n-1);
end
end
