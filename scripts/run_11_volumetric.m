%RUN_11_VOLUMETRIC  Does using the neighbouring B-scans as evidence help?
%
%   NEW FILE -- changes nothing that already existed.
%   Outputs: fig_22_volumetric.png, table_10_volumetric.csv
%
% THE QUESTION
%   A human grader does not look at one slice in isolation: a cyst is a
%   three-dimensional pocket, so if fluid is present in the slices either side
%   it is far more likely to be real.  Can we use that?
%
% WHAT THIS IS NOT
%   We already tested and rejected 3-D MEDIAN DENOISING (mixing the images of
%   k-1, k, k+1 before segmenting): Dice 0.426 -> 0.341, because smearing the
%   images across slices moves the mask away from the tracing drawn on slice
%   k.  See docs/FINDINGS.md section 7d.
%
%   This is the other way round: segment every slice INDEPENDENTLY, keep the
%   images separate, and use the neighbouring RESULTS as supporting evidence.
%
% WHAT WE CAN AND CANNOT VALIDATE
%   Only 110 of the 610 Duke B-scans were graded, and the graded ones are not
%   adjacent (typically every 5th slice).  So the neighbours k-1 and k+1
%   exist as IMAGES but carry no ground truth.  That is fine here: they are
%   used only as context, and scoring still happens on the graded slice.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

T = buildEvalSet(); E = T(T.evaluable,:);
n = height(E);
fprintf('Volumetric consistency test on %d evaluable B-scans\n', n);
fprintf('(each scan needs 3 segmentations, so this takes ~3x a normal run)\n\n');

% Each variant is a different rule for combining the three independent
% results.  'none' reproduces the current 2-D pipeline exactly and is the
% reference row.
variants = {
    'none',      struct('mode','none')
    'intersect', struct('mode','intersect','tolerancePx',3)
    'majority',  struct('mode','majority','tolerancePx',3)
    'union',     struct('mode','union','tolerancePx',3)
    'support 0.15', struct('mode','support','tolerancePx',3,'support',0.15)
    'support 0.30', struct('mode','support','tolerancePx',3,'support',0.30)
    'support 0.50', struct('mode','support','tolerancePx',3,'support',0.50)
    'support 0.70', struct('mode','support','tolerancePx',3,'support',0.70)
};

rows = {};
allDice = cell(size(variants,1),1);

for v = 1:size(variants,1)
    name = variants{v,1}; opt = variants{v,2};
    d = nan(n,1); pr = d; rc = d; tt = d;
    for i = 1:n
        t0 = tic;
        r  = segmentCME3D(E.subject(i), E.bscan(i), cfg, opt);
        tt(i) = toc(t0);
        S  = r.scan;
        valid = r.roi.mask & S.gradedMask;
        Ei = evalSegmentation(r.mask & valid, S.fluid1, valid);
        d(i) = Ei.dice; pr(i) = Ei.precision; rc(i) = Ei.recall;
    end
    allDice{v} = d;
    rows(end+1,:) = {name, mean(d,'omitnan'), mean(pr,'omitnan'), ...
                     mean(rc,'omitnan'), mean(tt)}; %#ok<SAGROW>
    fprintf('  %-14s Dice %.3f  Prec %.3f  Rec %.3f  (%.2f s/scan)\n', ...
        name, rows{end,2}, rows{end,3}, rows{end,4}, rows{end,5});
end

Tab = cell2table(rows, 'VariableNames', ...
    {'Rule','Dice','Precision','Recall','SecPerScan'});
writetable(Tab, fullfile(p.tables,'table_10_volumetric.csv'));

% ---- is any improvement real, or just noise? -------------------------
% 78 paired measurements, so a paired test on the per-scan Dice is the right
% check.  Without it, a 0.01 change means nothing.
base = allDice{1};
fprintf('\nPaired comparison against the current 2-D pipeline:\n');
sig = strings(size(variants,1),1); sig(1) = "reference";
for v = 2:size(variants,1)
    delta = allDice{v} - base;
    good  = isfinite(delta);
    if any(good)
        [~, pval] = ttest(delta(good));
        sig(v) = sprintf('%+.3f  (p = %.3f)', mean(delta(good)), pval);
    end
    fprintf('  %-14s %s\n', variants{v,1}, sig(v));
end
Tab.VsBaseline = sig;
writetable(Tab, fullfile(p.tables,'table_10_volumetric.csv'));

%% ----------------------------------------------------------------------
%  WHY does neighbour evidence not help?  The decisive test.
%  ----------------------------------------------------------------------
%  Neighbour agreement can only remove errors that are RANDOM from slice to
%  slice.  If our false positives are systematic -- the same dark structure
%  reappearing in every slice -- they will be just as well supported by the
%  neighbours as the true cysts, and no consistency rule can separate them.
%  So: measure the neighbour support of true and false components separately.
fprintf('\nDiagnostic: is neighbour support able to tell TP from FP at all?\n');
supTP = []; supFP = [];
for i = 1:n
    r = segmentCME3D(E.subject(i), E.bscan(i), cfg, ...
                     struct('mode','support','tolerancePx',3,'support',0));
    S = r.scan; valid = r.roi.mask & S.gradedMask;
    [L, nc] = bwlabel(r.mask2D & valid, 8);
    st = regionprops(L,'PixelIdxList','Area');
    for c = 1:nc
        idx = st(c).PixelIdxList;
        sup = nnz(r.votes(idx) > 0) / numel(idx);
        ov  = nnz(S.fluid1(idx)) / numel(idx);
        if ov > 0.30, supTP(end+1) = sup; else, supFP(end+1) = sup; end %#ok<SAGROW>
    end
end
na = min(numel(supTP), 4000);
aucSup = mean(randsample(supTP,na) > randsample(supFP,min(numel(supFP),na))');
fprintf('  true components : neighbour support %.3f (median)\n', median(supTP));
fprintf('  false components: neighbour support %.3f (median)\n', median(supFP));
fprintf('  separability (AUC) = %.3f   [0.5 = no discrimination at all]\n', mean(aucSup(:)));
fprintf('  n = %d true, %d false components\n', numel(supTP), numel(supFP));
save(fullfile(p.tables,'volumetricSupport.mat'), 'supTP','supFP');

%% ----------------------------------------------------------------------
f = figure('Position',[20 20 1500 640]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {'Using the neighbouring B-scans as evidence', ...
     'each slice is segmented independently; only the combination rule changes'});

nexttile;
b = bar(Tab.Dice,'FaceColor','flat');
b.CData = repmat([0.35 0.55 0.8], height(Tab), 1);
b.CData(1,:) = [0.45 0.45 0.45];                 % the 2-D reference
set(gca,'XTickLabel',Tab.Rule,'XTickLabelRotation',30);
ylabel('Dice'); ylim([0 0.6]); grid on
yline(Tab.Dice(1),'k--','LineWidth',1.5);
text(1:height(Tab), Tab.Dice, compose('%.3f',Tab.Dice), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
title('overlap with expert A');

nexttile;
% The decisive panel.  Neighbour agreement can only remove errors that vary
% from slice to slice.  These two distributions show it cannot tell our true
% components from our false ones, because BOTH are backed by the neighbours:
% the false positives are systematic structures that reappear in every slice.
edges = 0:0.05:1;
histogram(supTP, edges, 'Normalization','probability', ...
    'FaceColor',[0.20 0.60 0.30], 'FaceAlpha',0.65); hold on
histogram(supFP, edges, 'Normalization','probability', ...
    'FaceColor',[0.85 0.30 0.25], 'FaceAlpha',0.65);
grid on; xlabel('fraction of the component backed by a neighbouring slice');
ylabel('fraction of components');
legend({sprintf('TRUE cysts (median %.2f)', median(supTP)), ...
        sprintf('FALSE positives (median %.2f)', median(supFP))}, ...
    'Location','northwest');
title(sprintf(['both are ~95%% supported, so this cannot separate them' newline ...
    'separability AUC = %.3f   (0.5 = no discrimination)'], mean(aucSup(:))));

figStyle(f,'fig_22_volumetric');

fprintf('\nDone. Table: table_10_volumetric.csv,  figure: fig_22_volumetric.png\n');
