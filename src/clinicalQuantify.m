function C = clinicalQuantify(mask, roi, cfg)
%CLINICALQUANTIFY  Stage 5 -- OUR ADDITION.  Turn a contour into numbers a
%                  clinician can act on.
%
% Neither Liu et al. (2021) nor Lou et al. (2020) does this: both stop at the
% contour.  This module converts the segmentation into
%   * total intraretinal fluid area (px^2 and mm^2),
%   * cyst count,
%   * per-cyst and pooled equivalent-circle diameter (um),
%   * central subfield thickness (um), the quantity commercial OCT software
%     actually reports,
% and then hands those to gradeSeverity for a mild/moderate/severe label.
%
% INPUTS
%   mask  MxN logical final CME segmentation
%   roi   struct from extractROI (needed for retinal thickness / CST)
%   cfg   config struct (pixel pitches live here)
%
% OUTPUT (struct)
%   C.areaPx, C.areaMm2          total fluid area
%   C.cystCount                  number of connected components
%   C.perCyst                    table: areaPx, areaMm2, diamUm, widthUm,
%                                heightUm, centroidRow, centroidCol
%   C.meanDiamUm, C.maxDiamUm    pooled statistics
%   C.cstUm                      central subfield thickness
%   C.retinalThicknessUm         mean thickness over the whole B-scan
%   C.fluidFraction              fluid area / retinal area (unitless)
%   C.grade, C.gradeBasis        from gradeSeverity
%
% ANISOTROPIC PIXELS -- IMPORTANT
%   Duke/Spectralis pixels are NOT square: 3.87 um axially but ~11.46 um
%   laterally, a 3:1 aspect ratio.  So an "equivalent circle diameter"
%   computed in PIXELS would be meaningless.  We convert area to mm^2 FIRST
%   and derive the diameter from the physical area,
%        d = 2*sqrt(A/pi),
%   which is aspect-ratio correct.

px = cfg.pixelPitchAxialUm;      % um per row
py = cfg.pixelPitchLateralUm;    % um per column
pxAreaMm2 = cfg.pixelAreaMm2;

% ---------------------------------------------------------------------
% Totals
% ---------------------------------------------------------------------
C.areaPx   = nnz(mask);
C.areaMm2  = C.areaPx * pxAreaMm2;

% ---------------------------------------------------------------------
% Per-cyst measurements
% ---------------------------------------------------------------------
% bwconncomp finds the 8-connected components without building a full label
% image, and regionprops measures each one.  'MajorAxisLength' fits an
% ellipse with the same second moments as the region, which is the standard
% way to describe an elongated cyst.
CCs = bwconncomp(mask, 8);
S   = regionprops(CCs, 'Area','BoundingBox','Centroid','MajorAxisLength','MinorAxisLength');
C.cystCount = CCs.NumObjects;

n = C.cystCount;
areaPx = zeros(n,1); areaMm2 = zeros(n,1); diamUm = zeros(n,1);
widthUm = zeros(n,1); heightUm = zeros(n,1);
cRow = zeros(n,1); cCol = zeros(n,1);

for k = 1:n
    areaPx(k)  = S(k).Area;
    areaMm2(k) = areaPx(k) * pxAreaMm2;
    % equivalent-circle diameter from the PHYSICAL area (handles anisotropy)
    diamUm(k)  = 2 * sqrt((areaMm2(k)*1e6) / pi);      % mm^2 -> um^2 -> um
    widthUm(k)  = S(k).BoundingBox(3) * py;            % lateral extent
    heightUm(k) = S(k).BoundingBox(4) * px;            % axial extent
    cRow(k) = S(k).Centroid(2);
    cCol(k) = S(k).Centroid(1);
end

C.perCyst = table(areaPx, areaMm2, diamUm, widthUm, heightUm, cRow, cCol, ...
    'VariableNames', {'areaPx','areaMm2','diamUm','widthUm','heightUm','centroidRow','centroidCol'});

if n > 0
    C.meanDiamUm = mean(diamUm);
    C.maxDiamUm  = max(diamUm);
    C.maxCystAreaMm2 = max(areaMm2);
else
    C.meanDiamUm = 0; C.maxDiamUm = 0; C.maxCystAreaMm2 = 0;
end

% ---------------------------------------------------------------------
% Retinal thickness and central subfield thickness
% ---------------------------------------------------------------------
thicknessPxPerCol = roi.obm - roi.ilm;
C.retinalThicknessUm = mean(thicknessPxPerCol) * px;

% Central subfield = a cfg.clinical.centralMm wide strip centred on the frame.
% (The fovea is at the centre of the Duke 61-line scan protocol.)
N = numel(thicknessPxPerCol);
halfCols = round((cfg.clinical.centralMm * 1000 / py) / 2);
lo = max(1, round(N/2) - halfCols);
hi = min(N, round(N/2) + halfCols);
C.cstUm = mean(thicknessPxPerCol(lo:hi)) * px;
C.centralColumns = [lo hi];

retinaAreaPx = nnz(roi.mask);
C.fluidFraction = C.areaPx / max(retinaAreaPx, 1);
C.retinaAreaMm2 = retinaAreaPx * pxAreaMm2;

% ---------------------------------------------------------------------
% Severity grade
% ---------------------------------------------------------------------
[C.grade, C.gradeBasis, C.gradeDetail] = gradeSeverity(C, cfg);
end
