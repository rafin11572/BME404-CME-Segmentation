function W = waveMaps(I, cfg, analysisMask)
%WAVEMAPS  The wave potential-energy equation + potential-energy correction
%          equation, evaluated for ONE canonical operating direction.
%
% CANONICAL DIRECTION
%   Everything in here assumes the operator sweeps DOWN the rows, i.e. the
%   "wave" travels in the +row direction.  "Ahead / in front of" a pixel means
%   a LARGER row index; "behind" means a smaller one.  All other operating
%   directions are obtained by rotating the image (see omniWaveContours.m),
%   which is exactly the direction-adjustment function of Liu 2021 Eq. (5).
%
%   Index mapping to the papers: the papers write I(i,j) with i = x (column)
%   and j = y (depth), and their difference I(i,j+1)-I(i,j) steps one pixel
%   DEEPER.  In MATLAB that is I(r+1,c)-I(r,c).  So the papers' second index
%   is our ROW index.  Lou 2020 Sec. II-E confirms this ("the ILM and IRPE are
%   segmented along the direction of the cols").
%
% INPUTS
%   I             MxN double, the preprocessed working image on 0..255
%   cfg           config struct from cmeConfig
%   analysisMask  MxN logical (optional). Used only to compute the alpha
%                 threshold from the retina rather than from the vitreous.
%
% OUTPUT (struct, every field MxN, for both segmentation and figures)
%   W.phi_g   gravitational potential energy, Eq. (6)  -- in [0,1]
%   W.v       central velocity in the 3x3 Z template, Eq. (3)  -- in [-1,1]
%   W.vq      wave-front velocity in the 3x2 Q template, Eq. (4) -- in [-1,1]
%   W.sigma   speed regulation factor, Eq. (8)  -- in [0,1]
%   W.phi_v   kinetic energy = v .* vq .* sigma, Eq. (7)
%   W.wave    wave potential energy = phi_v + phi_g, Eq. (3)/(7)
%   W.alpha   the scalar threshold that was used
%   W.IT      logical "wave area" map, Lou 2020 Eq. (8)
%   W.K       geometric discriminant, Eq. (10)
%   W.C       gray discriminant,      Eq. (11)
%   W.Pen     logical "coastline" map = (K>1) & (C>1), Eq. (9)
%   W.bnd     logical boundary map = junction of coastline and wave areas
%
% WHY THIS FINDS CYST WALLS
%   Sweeping downwards, the operator fires where the image goes DARK -> BRIGHT
%   (C>1 demands the three pixels ahead to be brighter than the three behind).
%   A cyst is a dark lake inside brighter retina, so a downward sweep picks up
%   its FLOOR, an upward sweep its ROOF, and the two sideways sweeps its
%   walls.  That is exactly why four directions are needed to close the
%   contour, and it is the whole point of the 2021 paper.

D = double(I);
[M,N] = size(D);
if nargin < 3 || isempty(analysisMask), analysisMask = true(M,N); end
epsd = cfg.wave.epsDiv;

MAXg = max(D(:));
if MAXg <= 0, MAXg = 1; end

% =====================================================================
% phi_g  --  gravitational potential energy,  Liu 2021 Eq. (6)
%   phi_g = sum_Z I(i',j') * g(i',j')  /  MAX
% g is the Gaussian weighting function over the 3x3 Z template.  fspecial
% builds that kernel and normalises it to sum to 1, so the numerator is just a
% Gaussian-weighted local mean -- a noise-robust stand-in for the raw gray
% value ("gh" in the fluid equation, with g=1 and h=gray).  Dividing by the
% image maximum puts phi_g in [0,1] as the paper requires.
% =====================================================================
gk = fspecial('gaussian', [3 3], cfg.wave.gaussSigma);
W.phi_g = imfilter(D, gk, 'replicate', 'same') / MAXg;

% =====================================================================
% Forward differences along the sweep direction.
%   d(r,c) = I(r+1,c) - I(r,c)
% This is the difference approximation to the first derivative of displacement
% that the papers use in place of the fluid velocity (Lou 2020 Sec. II-C).
% The last row has no successor, so it is padded with zero.
% =====================================================================
d = [diff(D,1,1); zeros(1,N)];

% ---- v : normalised signed difference inside the 3x3 Z template, Eq.(3) ----
% A 3x3 Z template centred on (r,c) fully contains the differences
% d(r-1,c') and d(r,c') for c' = c-1..c+1, i.e. six signed differences.
% Correlation with this 3x3 mask of ones sums exactly those six.
kZ = [1 1 1; 1 1 1; 0 0 0];
sumZ    = imfilter(d,      kZ, 'replicate', 'corr');
sumAbsZ = imfilter(abs(d), kZ, 'replicate', 'corr');
W.v = sumZ ./ max(sumAbsZ, epsd);      % in [-1,1]; ->1 when the ramp is
                                       % monotonic and steep in one direction

% ---- vq : same thing inside the 3x2 Q template placed AHEAD of Z, Eq.(4) ---
% Q is 3 wide (columns c-1..c+1) and 2 deep, sitting just in front of the Z
% template, so it spans rows r+2 and r+3.  The single difference it contains
% per column is d(r+2,c').  A 5x3 correlation kernel with ones on its last row
% reaches exactly r+2 when centred at r.
kQ = zeros(5,3); kQ(5,:) = 1;
sumQ    = imfilter(d,      kQ, 'replicate', 'corr');
sumAbsQ = imfilter(abs(d), kQ, 'replicate', 'corr');
W.vq = sumQ ./ max(sumAbsQ, epsd);
% vq is the "wind vane" of v: kinetic energy v*vq only stays large when the
% ramp continues ahead of the template too.  An isolated speckle spike ("reef")
% gives a large v but a random vq, so the product collapses.  This is how the
% method kills speckle without ever blurring it away.

% =====================================================================
% sigma  --  speed regulation factor, Eq. (8)
%   sigma = exp( -( I / (255*delta) )^2 ),   delta = median(Z)/mean(Z)
% medfilt2 gives the median of the 3x3 Z template at every pixel; correlation
% with ones(3)/9 gives its mean.  In a flat region median~mean so delta~1 and
% sigma = exp(-(I/255)^2): bright pixels get sigma -> exp(-1) = 0.37, dark
% pixels sigma -> 1.  That damps the kinetic energy inside the hyper-reflective
% retinal layers, which would otherwise be over-segmented.
% =====================================================================
md   = medfilt2(D, [3 3], 'symmetric');
mn   = imfilter(D, ones(3)/9, 'replicate', 'corr');
delta = md ./ max(mn, epsd);
W.sigma = exp( -( D ./ (255 * max(delta, epsd)) ).^2 );

% ---- kinetic energy and total wave potential energy, Eq. (7) and (3) ----
W.phi_v = W.v .* W.vq .* W.sigma;
W.wave  = W.phi_v + W.phi_g;

% =====================================================================
% Binarise the wave map into the "wave area" IT.
% =====================================================================
[W.IT, W.alpha] = thresholdWaveMap(W.wave, D, analysisMask, cfg);

% =====================================================================
% Potential energy correction equation, Lou 2020 Eqs. (9)-(11)
%   K = Kb/Kc = ( I(j-1) - I(j-3) ) / ( I(j+2) - I(j) )     geometric
%   C = Sc/Sb = ( I(j)+I(j+1)+I(j+2) ) / ( I(j-1)+I(j-2)+I(j-3) )  gray
%   Penergy = 1  iff  K > 1  AND  C > 1
% The offsets really are asymmetric in the paper (behind uses -1 and -3, ahead
% uses 0 and +2) and we implement them literally.
% =====================================================================
Im3 = shiftRows(D,-3); Im2 = shiftRows(D,-2); Im1 = shiftRows(D,-1);
Ip1 = shiftRows(D, 1); Ip2 = shiftRows(D, 2);

Kb = Im1 - Im3;          % slope of the three points BEHIND the measuring point
Kc = Ip2 - D;            % slope of the measuring point and two points AHEAD

switch lower(cfg.wave.correctionMode)
    case 'ratio'
        % Literal printed form.  Where Kc is ~0 the ratio is +-Inf, which is
        % the correct limit (an infinitely steeper slope behind than ahead).
        W.K = Kb ./ Kc;
        W.K(abs(Kc) < epsd) = sign(Kb(abs(Kc) < epsd)) .* Inf;
        Kok = W.K > 1;
    case 'slope'
        % Sign-safe form, matching the paper's prose "if K is greater than 1
        % (Kb > Kc)".  These disagree whenever Kc < 0.
        W.K = Kb - Kc;
        Kok = Kb > Kc;
    otherwise
        error('waveMaps:correctionMode','Unknown correctionMode');
end

Sc = D + Ip1 + Ip2;      % sum of the measuring point and its two front points
Sb = Im1 + Im2 + Im3;    % sum of the three points behind it
W.C = Sc ./ max(Sb, epsd);
Cok = W.C > 1;

W.Pen = Kok & Cok;       % the "coastline" area

% =====================================================================
% The real boundary sits at the JUNCTION of the coastline and the wave area
% (Lou 2020 Sec. II-D): "the retinal boundary line can be obtained by
% extracting the points at the junction of the coastline area and the wave
% area".
%
% What a junction IS, concretely.  Sweep downwards into the floor of a cyst:
%     inside the cyst   -- dark, so phi_g is low  -> wave low  -> IT = 0
%     crossing the wall -- gray ramps up, v and vq both -> 1, phi_v peaks,
%                          phi_g rising            -> wave high -> IT = 1
%     inside the tissue -- bright but flat, phi_v ~ 0, phi_g high -> IT = 1
% So the wall is the pixel at which IT SWITCHES ON.  Taking every pixel of the
% wave area (or every coastline pixel) would mark ~60% of the retina; taking
% the leading edge marks a thin line, which is what Fig. 5 of the 2021 paper
% actually shows.  The correction equation Pen then vets that edge, rejecting
% the ones whose local slope pattern is not a true boundary.
%
%   'edge' (default) : Pen AND the leading edge of the wave area
%   'all'            : Pen anywhere that touches the wave area -- the naive
%                      literal reading; kept so the difference can be shown.
if ~isfield(cfg.wave,'boundaryMode'), cfg.wave.boundaryMode = 'edge'; end
switch lower(cfg.wave.boundaryMode)
    case 'edge'
        leading = W.IT & ~shiftRows(W.IT,-1);   % IT on here, off one step back
        W.bnd = W.Pen & leading;
    case 'all'
        W.bnd = W.Pen & (W.IT | shiftRows(W.IT,-1));
    otherwise
        error('waveMaps:boundaryMode','Unknown boundaryMode');
end
W.bnd = W.bnd & analysisMask;
end


% =========================================================================
function S = shiftRows(X, k)
%SHIFTROWS  S(r,c) = X(r+k,c), with the image replicated at the borders.
%   k > 0 looks AHEAD (deeper), k < 0 looks BEHIND (shallower).
[M,~] = size(X);
idx = min(max((1:M)' + k, 1), M);   % clamp instead of wrapping: circshift
S = X(idx, :);                      % would make the top of the image "ahead"
end                                 % of the bottom, which is physically wrong


% =========================================================================
function [IT, alpha] = thresholdWaveMap(waveMap, D, mask, cfg)
%THRESHOLDWAVEMAP  Turn the continuous wave map into the binary wave area.
%
%   'alpha'   [PAPER] Lou 2020 Eq.(8):  IT = wave > alpha,  alpha = meanB/meanT
%                     with meanB the mean gray of the background area and
%                     meanT the mean gray of the target area.
%   'otsu'    [OURS]  global Otsu on the wave map
%   'niblack' [OURS]  Niblack local threshold, T = localmean + k*localstd

switch lower(cfg.wave.threshold)

    case 'alpha'
        % Split the analysed region into "background" and "target" by Otsu on
        % the GRAY image (not the wave map) -- that is what meanB / meanT mean
        % in the paper: dark background vs bright tissue.
        g = D(mask);
        if isempty(g), g = D(:); end
        lvl = graythresh(uint8(g)) * 255;   % graythresh returns a normalised
                                            % level; see niblackThreshold.m for
                                            % what Otsu actually optimises
        meanB = mean(g(g <= lvl));
        meanT = mean(g(g >  lvl));
        if isnan(meanB) || isnan(meanT) || meanT == 0
            alpha = 0.5;
        else
            alpha = meanB / meanT;
        end
        IT = waveMap > alpha;

    case 'otsu'
        % Normalise the wave map to [0,1] then let Otsu choose the level.
        % graythresh implements Otsu's method: it searches every possible
        % threshold and keeps the one that MAXIMISES the between-class variance
        % of the two resulting intensity groups (equivalently, minimises the
        % within-class variance).  It is fully automatic and needs no tuning,
        % which is the point of our modification -- alpha above is a fixed
        % global ratio that does not adapt to an unusually bright or dark scan.
        wn = normalise01(waveMap);
        alpha = graythresh(wn);
        IT = wn > alpha;

    case 'niblack'
        wn = normalise01(waveMap);
        [IT, alpha] = niblackThreshold(wn, cfg.wave.niblackWindow, cfg.wave.niblackK);

    otherwise
        error('waveMaps:threshold','Unknown wave threshold "%s"', cfg.wave.threshold);
end
end


% =========================================================================
function y = normalise01(x)
lo = min(x(:)); hi = max(x(:));
if hi <= lo, y = zeros(size(x)); else, y = (x - lo) / (hi - lo); end
end
