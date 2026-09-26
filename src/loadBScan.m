function S = loadBScan(subjectId, k, halfWidth)
%LOADBSCAN  Load ONE B-scan and all of its annotations from the Duke DME set.
%
% INPUTS
%   subjectId  integer 1..10  (Subject_01.mat ... Subject_10.mat)
%   k          integer 1..61  index along the third axis (the B-scan number)
%   halfWidth  how many NEIGHBOURING B-scans to load on each side (default 1).
%              Speckle is uncorrelated between B-scans but anatomy is not, so a
%              median across k-1,k,k+1 suppresses speckle without blurring the
%              cyst walls the way an in-plane filter does.  See S.stack.
%
% OUTPUT  (struct)
%   S.img        496x768 double, gray values 0..255 (the raw B-scan)
%   S.fluid1     496x768 logical, expert A manual fluid mask  (manualFluid1)
%   S.fluid2     496x768 logical, expert B manual fluid mask  (manualFluid2)
%   S.autoFluid  496x768 logical, the dataset's own automatic result
%                (NOT ground truth -- only ever used as an extra comparator)
%   S.layers1    8x768 double, expert A boundary traces (NaN where untraced)
%   S.layers2    8x768 double, expert B boundary traces
%   S.subjectId, S.bscan, S.hasFluid, S.hasLayers
%
% The 8 layer traces are, in order:
%   1 ILM, 2 NFL/GCL, 3 IPL/INL, 4 INL/OPL, 5 OPL/ONL, 6 ISM/ISE, 7 OS/RPE, 8 BM
% They are MATLAB 1-based ROW indices, one value per A-scan (column).
%
% IMPLEMENTATION NOTE
%   We use matfile() rather than load().  matfile gives *partial* access to a
%   variable inside a .mat file, so pulling one 496x768 slice out of a
%   496x768x61 array reads only that slice off disk instead of all 20 MB.
%   That is what makes the GUI feel instant when you scrub through B-scans.

if nargin < 3 || isempty(halfWidth), halfWidth = 1; end
p = dukePaths();
% Prefer the v7.3 copy if it exists: only that format supports real partial
% loading, so reading one slice costs one slice instead of the whole 180 MB
% variable.  Run convertToV73 once to create it.  Falls back transparently.
f = fullfile(p.data, 'v73', sprintf('Subject_%02d.mat', subjectId));
if exist(f,'file') ~= 2
    f = fullfile(p.data, sprintf('Subject_%02d.mat', subjectId));
end
assert(exist(f,'file')==2, 'loadBScan:missingFile', ...
    'Cannot find Subject_%02d.mat in %s', subjectId, p.data);

m = matfile(f);                       % lazy handle, nothing read yet

S.subjectId = subjectId;
S.bscan     = k;
S.img       = double(m.images(:,:,k));

% The manual masks are stored as double with NaN outside the graded window and
% 1 inside a cyst.  ">0" turns that into a clean logical mask and also kills
% the NaNs (NaN>0 is false).
S.fluid1    = m.manualFluid1(:,:,k) > 0;
S.fluid2    = m.manualFluid2(:,:,k) > 0;
S.autoFluid = m.automaticFluidDME(:,:,k) > 0;

S.layers1   = squeeze(m.manualLayers1(:,:,k));   % 8 x 768
S.layers2   = squeeze(m.manualLayers2(:,:,k));

% Neighbouring B-scans, centre slice at index S.centre.
ks = max(1,k-halfWidth):min(size(m,'images',3), k+halfWidth);
S.stack  = double(m.images(:,:,ks));
S.centre = find(ks == k);

% Was this B-scan graded at all?  The Duke experts annotated only 110 of the
% 610 B-scans; on the other 500 manualFluid1 is entirely NaN.  Those scans are
% UNGRADED, not fluid-free, so "we segmented something and the expert did not"
% means nothing there.  The layer traces mark which scans were graded.
S.isGraded  = any(isfinite(S.layers1(:)));
S.hasFluid  = any(S.fluid1(:));
S.hasLayers = any(isfinite(S.layers1(:)));

% The Duke graders did not annotate the full 768-A-scan width: they worked
% inside a central window, and the layer traces mark its extent.  We verified
% across all 78 fluid-positive B-scans that every annotated fluid pixel lies
% strictly inside that window, so anything we segment outside it cannot be
% checked against a grader and must not be counted as a false positive.
% S.gradedCols is that window, and the evaluation is restricted to it.
S.gradedCols = any(isfinite(S.layers1), 1);      % 1x768 logical
if ~any(S.gradedCols), S.gradedCols = true(1, size(S.img,2)); end
S.gradedMask = repmat(S.gradedCols, size(S.img,1), 1);
end
