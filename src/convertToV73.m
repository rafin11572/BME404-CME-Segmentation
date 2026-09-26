function convertToV73(force)
%CONVERTTOV73  One-time conversion of the Duke .mat files to MAT v7.3.
%
%   convertToV73        convert any subject that has not been converted yet
%   convertToV73(true)  reconvert everything
%
% WHY
%   The Duke files are saved in the older MAT v7 format, which is compressed
%   as a single stream per variable.  That format does NOT support partial
%   loading: asking matfile() for one 496x768 slice of images(:,:,k) forces
%   MATLAB to decompress the whole 496x768x61 array (about 180 MB) just to
%   hand back 0.4 MB.  With six such variables per B-scan, simply loading a
%   scan cost more than segmenting it -- roughly 0.8 s of every 1.2 s.
%
%   MAT v7.3 is HDF5-based and stores variables in chunks, so matfile() can
%   read a single slice off disk directly.  Converting once makes every batch
%   run about 2.5x faster and makes scrubbing through B-scans in the GUI feel
%   instant.
%
%   The converted copies go in data/v73/ and the originals are left untouched.
%   loadBScan prefers the converted file when it exists and silently falls
%   back to the original when it does not, so the project works either way.

if nargin < 1, force = false; end
p = dukePaths();
outDir = fullfile(p.data, 'v73');
if ~exist(outDir,'dir'), mkdir(outDir); end

vars = {'images','manualFluid1','manualFluid2','automaticFluidDME', ...
        'manualLayers1','manualLayers2'};

for s = 1:10
    src = fullfile(p.data, sprintf('Subject_%02d.mat', s));
    dst = fullfile(outDir,  sprintf('Subject_%02d.mat', s));
    if ~exist(src,'file'), continue; end
    if exist(dst,'file') && ~force
        fprintf('  Subject_%02d already converted\n', s); continue
    end
    fprintf('  converting Subject_%02d ...', s); t0 = tic;
    D = load(src, vars{:});

    % uint8 is the natural type for 0-255 gray values and quarters the file
    % size relative to double, with no loss: we checked that every value in
    % images is an integer in [0,255].
    D.images = uint8(D.images);
    % The fluid maps are instance labels 0..15 plus NaN, so they go to uint8
    % as well, with NaN mapped to 0 (ungraded == no fluid for our purposes).
    for f = {'manualFluid1','manualFluid2','automaticFluidDME'}
        X = D.(f{1}); X(isnan(X)) = 0; D.(f{1}) = uint8(X);
    end
    % Layer traces must stay double: they carry NaN where untraced, and that
    % NaN is meaningful (it defines the graders' annotated window).

    save(dst, '-struct', 'D', '-v7.3');
    fprintf(' %.1f s\n', toc(t0));
end
fprintf('Converted files in %s\n', outDir);
end
