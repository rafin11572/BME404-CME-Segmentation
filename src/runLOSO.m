function R = runLOSO(baseCfg, grid, verbose)
%RUNLOSO  Leave-one-subject-out validation of the tuned parameters.
%
%   R = runLOSO(cmeConfig('modified'), grid)
%
% WHY THIS EXISTS
%   Every threshold in this pipeline was chosen by looking at how well it
%   scored on the Duke B-scans.  Reporting that same score as the pipeline's
%   accuracy would be circular: it measures how well we fitted the test set,
%   not how well the method generalises.
%
%   Leave-one-subject-out fixes that.  For each subject in turn we re-tune the
%   free parameters on the OTHER NINE subjects only, then score the held-out
%   subject with those parameters.  No B-scan is ever scored with parameters
%   that were chosen using that same patient's data.  This is the protocol
%   Rashno et al. (2018) use on this dataset, and it is the number that should
%   be quoted in the report.
%
% INPUTS
%   baseCfg  the configuration to start from (cmeConfig('modified'))
%   grid     struct array with fields .name (dotted config path) and .values
%            (cell array) -- the free parameters to re-tune on each fold
%   verbose  print per-fold progress (default true)
%
% OUTPUT (struct)
%   R.perScan     table of held-out scores for every evaluable B-scan
%   R.perFold     table: subject, the parameters chosen on the other nine,
%                 and the held-out Dice
%   R.summary     mean held-out Dice / precision / recall
%   R.fittedAll   the same metrics when tuned AND tested on everything, for
%                 comparison -- the gap between the two is the optimism

if nargin < 3, verbose = true; end

T = buildEvalSet(); E = T(T.evaluable,:);
subjects = unique(E.subject)';

% ---- cache every scan once -------------------------------------------
S = cell(height(E),1);
for i = 1:height(E)
    S{i} = loadBScan(E.subject(i), E.bscan(i), 1);
end

combos = buildCombos(grid);
if verbose
    fprintf('LOSO: %d subjects, %d parameter combinations, %d B-scans\n', ...
        numel(subjects), numel(combos), height(E));
end

% ---- score every combination on every scan ONCE ----------------------
% The folds differ only in which scans are averaged, so the expensive part
% (segmenting) is done once per combination and reused by every fold.
D = nan(numel(combos), height(E));
P = D; Rc = D;
for c = 1:numel(combos)
    cfg = applyCombo(baseCfg, grid, combos{c});
    for i = 1:height(E)
        res = segmentCME(S{i}.img, cfg, S{i}.stack);
        v   = res.roi.mask & S{i}.gradedMask;
        Ei  = evalSegmentation(res.mask & v, S{i}.fluid1, v);
        D(c,i) = Ei.dice; P(c,i) = Ei.precision; Rc(c,i) = Ei.recall;
    end
    if verbose, fprintf('  combo %d/%d done (mean Dice %.3f)\n', c, numel(combos), mean(D(c,:),'omitnan')); end
end

% ---- sanity check: did every grid dimension actually DO anything? -----
% A parameter that the current configuration never reads (for example an
% ROI parameter that only applies to a band mode we are not using) silently
% contributes nothing: the grid looks larger than the search really was, and
% the "chosen parameters" string ends up quoting a meaningless tie-break.
% This caught exactly that mistake once, so it is checked automatically.
R.inertParameters = strings(0,1);
for g = 1:numel(grid)
    others = setdiff(1:numel(grid), g);
    key = zeros(numel(combos), max(numel(others),1));
    for c = 1:numel(combos)
        if isempty(others), key(c,1) = 1; else, key(c,:) = combos{c}(others); end
    end
    [~,~,grp] = unique(key, 'rows');
    inert = true;
    for q = unique(grp)'
        rows = D(grp == q, :);
        if size(rows,1) > 1 && any(max(rows,[],1) - min(rows,[],1) > 1e-12)
            inert = false; break
        end
    end
    if inert
        R.inertParameters(end+1) = string(grid(g).name);
        warning('runLOSO:inertParameter', ...
            ['"%s" had NO effect on any result -- the current configuration ' ...
             'never reads it. Remove it from the grid or it will inflate the ' ...
             'apparent size of the search.'], grid(g).name);
    end
end

% ---- fold loop --------------------------------------------------------
heldDice = nan(height(E),1); heldPrec = heldDice; heldRec = heldDice;
foldRows = {};
for s = subjects
    test  = (E.subject == s);
    train = ~test;
    % choose the combination that is best on the OTHER subjects
    trainScore = mean(D(:,train), 2, 'omitnan');
    [~, bestC] = max(trainScore);
    heldDice(test) = D(bestC, test);
    heldPrec(test) = P(bestC, test);
    heldRec(test)  = Rc(bestC, test);
    foldRows(end+1,:) = {s, comboLabel(grid, combos{bestC}), ...
        trainScore(bestC), mean(D(bestC,test),'omitnan')}; %#ok<AGROW>
end

R.perScan = table(E.subject, E.bscan, heldDice, heldPrec, heldRec, ...
    'VariableNames', {'subject','bscan','dice','precision','recall'});
R.perFold = cell2table(foldRows, ...
    'VariableNames', {'heldOutSubject','parametersChosen','trainDice','heldOutDice'});
R.summary = table(mean(heldDice,'omitnan'), mean(heldPrec,'omitnan'), mean(heldRec,'omitnan'), ...
    'VariableNames', {'Dice','Precision','Recall'});

% the optimistic number: best single combination over everything
allScore = mean(D, 2, 'omitnan');
[bestAll, iAll] = max(allScore);
R.fittedAll = table(bestAll, mean(P(iAll,:),'omitnan'), mean(Rc(iAll,:),'omitnan'), ...
    'VariableNames', {'Dice','Precision','Recall'});
R.fittedAllParams = comboLabel(grid, combos{iAll});
R.optimism = bestAll - mean(heldDice,'omitnan');

% ---- spread across folds ---------------------------------------------
% Ten folds over ten patients is a small sample, so the mean on its own is
% not an honest summary.  Report the spread with it.
fd = R.perFold.heldOutDice;
R.foldSpread = table(mean(fd), std(fd), min(fd), max(fd), ...
    'VariableNames', {'meanFoldDice','stdFoldDice','minFoldDice','maxFoldDice'});
R.effectiveCombos = size(unique(D,'rows'), 1);

% =====================================================================
% WHAT THIS NUMBER DOES AND DOES NOT MEASURE  -- read before quoting it
% =====================================================================
% It measures the optimism of the FINAL PARAMETER SELECTION step only: for
% each fold the values are picked using nine subjects and scored on the
% tenth, so no patient's own data chose the thresholds used on them.
%
% It does NOT measure the optimism of the whole design process.  The pipeline
% structure, the choice of which parameters to expose, and the ranges swept
% were all decided while looking at the entire dataset.  A fully honest
% estimate would need a locked-away test set that was never examined, which
% this dataset (10 patients, 78 usable B-scans) is too small to provide.
% So treat the held-out figure as an upper bound on how well the method would
% do on a new patient from the same scanner, not as a guarantee.
R.caveat = ['Measures the optimism of the final parameter selection only. ' ...
            'The pipeline structure and the swept ranges were chosen using ' ...
            'the whole dataset, so the true optimism is larger than reported.'];
end


% =========================================================================
function combos = buildCombos(grid)
n = arrayfun(@(g) numel(g.values), grid);
total = prod(n);
combos = cell(1,total);
for t = 1:total
    idx = cell(1,numel(grid));
    [idx{:}] = ind2sub(n, t);
    combos{t} = cell2mat(idx);
end
end

function cfg = applyCombo(cfg, grid, combo)
for g = 1:numel(grid)
    parts = strsplit(grid(g).name, '.');
    S = struct('type', repmat({'.'},1,numel(parts)), 'subs', parts);
    cfg = subsasgn(cfg, S, grid(g).values{combo(g)});
end
end

function s = comboLabel(grid, combo)
parts = cell(1,numel(grid));
for g = 1:numel(grid)
    v = grid(g).values{combo(g)};
    if ischar(v) || isstring(v), vs = char(v); else, vs = mat2str(v); end
    parts{g} = sprintf('%s=%s', grid(g).name, vs);
end
s = string(strjoin(parts, ', '));
end
