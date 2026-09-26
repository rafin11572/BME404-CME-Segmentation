function R = runBatch(cfg, scans, opts)
%RUNBATCH  Run the pipeline over a list of B-scans and score every one.
%
%   R = runBatch(cfg)                 all 78 fluid-positive B-scans
%   R = runBatch(cfg, scans)          scans = Kx2 [subject bscan]
%   R = runBatch(cfg, scans, opts)    opts.verbose, opts.keepMasks
%
% OUTPUT (struct)
%   R.perScan   table, one row per B-scan, with the metrics against expert A
%               and expert B and the processing time
%   R.summary   table in the shape of Table 1 of Liu et al. (2021):
%               rows = metric, columns = Expert A / Expert B / Average
%   R.cfg       the configuration used
%   R.masks     cell array of the output masks, if opts.keepMasks
%
% SCORING PROTOCOL (stated explicitly because it is easy to get wrong)
%   * Metrics are computed only on B-scans where expert A marked fluid --
%     Dice is undefined on an all-zero ground truth.  Scans with no fluid are
%     handled separately by runFalsePositiveCheck.
%   * TN is counted inside the retina (the ROI) only.  Over the whole
%     496x768 frame, accuracy would exceed 99% purely because most of the
%     image is empty vitreous.
%   * Pixels outside the graders' annotated A-scan window are excluded, since
%     no ground truth was ever drawn there.  See loadBScan.

if nargin < 2 || isempty(scans), scans = defaultEvalScans(); end
if nargin < 3, opts = struct(); end
if ~isfield(opts,'verbose'),   opts.verbose   = true;  end
if ~isfield(opts,'keepMasks'), opts.keepMasks = false; end

n = size(scans,1);
vars = {'subject','bscan','diceA','precA','recA','f1A','accA','jacA', ...
                          'diceB','precB','recB','f1B','accB','jacB', ...
        'areaPx','areaMm2','cystCount','meanDiamUm','cstUm','grade','timeSec'};
P = cell(n, numel(vars));
masks = cell(n,1);

for i = 1:n
    S = loadBScan(scans(i,1), scans(i,2));
    res = segmentCME(S.img, cfg, S.stack);

    valid = res.roi.mask;
    if cfg.eval.restrictToGraded, valid = valid & S.gradedMask; end
    if ~cfg.eval.restrictToROI,   valid = true(size(S.img)); end

    pred = res.mask & valid;
    EA = evalSegmentation(pred, S.fluid1, valid);
    EB = evalSegmentation(pred, S.fluid2, valid);

    C = res.clinical;
    P(i,:) = {scans(i,1), scans(i,2), ...
        EA.dice, EA.precision, EA.recall, EA.f1, EA.accuracy, EA.jaccard, ...
        EB.dice, EB.precision, EB.recall, EB.f1, EB.accuracy, EB.jaccard, ...
        C.areaPx, C.areaMm2, C.cystCount, C.meanDiamUm, C.cstUm, ...
        string(C.grade), res.timing.total};

    if opts.keepMasks, masks{i} = pred; end
    if opts.verbose && mod(i,10)==0
        fprintf('    %3d/%d scans\n', i, n);
    end
end

R.perScan = cell2table(P, 'VariableNames', vars);
R.summary = summarise(R.perScan);
R.cfg     = cfg;
if opts.keepMasks, R.masks = masks; R.scans = scans; end
end


% =========================================================================
function T = summarise(ps)
%SUMMARISE  Build the paper's Table 1 layout: metric x {Expert A, B, Average}.
%
% Both averaging conventions are reported.  The 2021 paper quotes Dice 0.811
% and F1 0.813 for what is algebraically the same quantity; the difference
% comes from averaging per image versus pooling all pixels, so we state which
% is which instead of presenting them as two different measures.
metric = {'Dice';'Precision';'Recall';'F1-Score';'Accuracy';'Jaccard'};
A = [mean(ps.diceA,'omitnan'); mean(ps.precA,'omitnan'); mean(ps.recA,'omitnan');
     mean(ps.f1A,'omitnan');   mean(ps.accA,'omitnan');  mean(ps.jacA,'omitnan')];
B = [mean(ps.diceB,'omitnan'); mean(ps.precB,'omitnan'); mean(ps.recB,'omitnan');
     mean(ps.f1B,'omitnan');   mean(ps.accB,'omitnan');  mean(ps.jacB,'omitnan')];
T = table(metric, A, B, (A+B)/2, 'VariableNames', {'Metric','ExpertA','ExpertB','Average'});
end


% =========================================================================
function scans = defaultEvalScans()
%DEFAULTEVALSCANS  The 78 B-scans on which expert A marked fluid.
T = buildEvalSet();
E = T(T.evaluable, :);
scans = [E.subject, E.bscan];
end
