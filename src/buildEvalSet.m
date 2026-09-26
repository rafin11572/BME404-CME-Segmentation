function T = buildEvalSet(forceRebuild)
%BUILDEVALSET  Index every B-scan in the Duke DME set and cache the result.
%
%   T = buildEvalSet()       load from cache if it exists, else build it
%   T = buildEvalSet(true)   force a rebuild
%
% OUTPUT
%   T  a table with one row per (subject, B-scan), columns:
%        subject, bscan, nFluid1, nFluid2, nAuto, hasLayers,
%        hasFluid1, hasFluid2, evaluable
%      "evaluable" means expert A marked at least one fluid pixel -- these are
%      the scans the segmentation metrics are computed on, exactly as the
%      papers do (you cannot define Dice on an all-zero ground truth).
%
% The scans WITHOUT fluid are kept in the table on purpose: they give us a
% false-positive-rate-on-healthy-scans robustness figure, which neither paper
% reports.
%
% Building the index reads all ten 20 MB files once (~1-2 min), so the result
% is cached in results/tables/evalIndex.mat.

if nargin < 1, forceRebuild = false; end
p = dukePaths();
cacheFile = fullfile(p.tables, 'evalIndex.mat');

if ~forceRebuild && exist(cacheFile,'file')
    S = load(cacheFile); T = S.T; return
end

rows = [];
for s = 1:10
    f = fullfile(p.data, sprintf('Subject_%02d.mat', s));
    if ~exist(f,'file')
        warning('buildEvalSet:missing','Skipping missing %s', f); continue
    end
    fprintf('  indexing Subject_%02d ...\n', s);
    D = load(f, 'manualFluid1','manualFluid2','automaticFluidDME','manualLayers1');

    % squeeze(sum(sum(X,1),2)) collapses each 496x768 slice to a single count,
    % giving a 61x1 vector of "how many fluid pixels in this B-scan".
    n1 = squeeze(sum(sum(D.manualFluid1 > 0, 1), 2));
    n2 = squeeze(sum(sum(D.manualFluid2 > 0, 1), 2));
    na = squeeze(sum(sum(D.automaticFluidDME > 0, 1), 2));
    nL = squeeze(sum(sum(isfinite(D.manualLayers1), 1), 2));

    nB = numel(n1);
    rows = [rows; table(repmat(s,nB,1), (1:nB)', n1, n2, na, nL>0, ...
        'VariableNames', {'subject','bscan','nFluid1','nFluid2','nAuto','hasLayers'})]; %#ok<AGROW>
end

rows.hasFluid1 = rows.nFluid1 > 0;
rows.hasFluid2 = rows.nFluid2 > 0;
rows.evaluable = rows.hasFluid1;
T = rows;

save(cacheFile, 'T');
writetable(T, fullfile(p.tables, 'evalIndex.csv'));

fprintf(['\nDuke DME index built:\n' ...
         '  %d B-scans total\n' ...
         '  %d with expert-A fluid (the evaluable set)\n' ...
         '  %d with expert-B fluid\n' ...
         '  %d with finite layer traces\n'], ...
        height(T), sum(T.hasFluid1), sum(T.hasFluid2), sum(T.hasLayers));
end
