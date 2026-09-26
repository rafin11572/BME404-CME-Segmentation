%RUN_05_FIGURES_OURS  The figures that show what WE added, none of which
%                     appear in either reference paper.
%
%   fig_09_baseline_vs_modified   side-by-side segmentation, same B-scans
%   fig_10_metrics_comparison     grouped bars: paper / our baseline / ours
%   fig_11_per_scan_dice          Dice for every one of the 78 evaluable scans
%   fig_12_interobserver          our agreement vs the experts' agreement
%   fig_13_failure_gallery        honest look at where the method breaks
%   fig_14_preprocessing_effect   why the CLAHE tile size matters

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

T = buildEvalSet(); E = T(T.evaluable,:);
scans = [E.subject, E.bscan];
L = load(fullfile(p.tables,'batchResults.mat'));  R = L.R;

cfgB = cmeConfig('baseline');
cfgM = cmeConfig('modified');

%% ======================================================================
%  fig_09  Baseline vs modified, side by side  (our clearest contribution)
%  ======================================================================
pick = [1 31; 5 24; 6 33; 8 19];
f = figure('Position',[10 10 1700 1120]);
tl = tiledlayout(f, size(pick,1), 3, 'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Published pipeline (Liu et al. 2021) vs our modified pipeline -- green = expert A, red = algorithm');

for i = 1:size(pick,1)
    S = loadBScan(pick(i,1), pick(i,2));
    rb = segmentCME(S.img, cfgB, S.stack);
    rm = segmentCME(S.img, cfgM, S.stack);
    valid = rm.roi.mask & S.gradedMask;
    Eb = evalSegmentation(rb.mask & valid, S.fluid1, valid);
    Em = evalSegmentation(rm.mask & valid, S.fluid1, valid);

    a0 = max(1, round(min(rm.roi.ilm))-25);
    a1 = min(size(S.img,1), round(max(rm.roi.obm))+25);
    cc = @(X) X(a0:a1,:,:);

    nexttile; imshow(vizOverlay(cc(S.img), {cc(S.fluid1)}, {[0 1 0]}, struct('lineWidth',2)));
    title(sprintf('S%02d k%02d  expert A', pick(i,1), pick(i,2)));

    nexttile; imshow(vizOverlay(cc(S.img), {cc(S.fluid1), cc(rb.mask & valid)}, ...
        {[0 1 0],[1 0 0]}, struct('fillAlpha',0.18,'lineWidth',2)));
    if nnz(rb.mask & valid) == 0
        text(15, 20, 'baseline output: EMPTY', 'Color',[1 0.5 0.3], 'FontSize',11, ...
            'FontWeight','bold','VerticalAlignment','top','BackgroundColor','k','Margin',3);
    end
    title(sprintf('baseline (paper-literal)   Dice %.2f', Eb.dice));

    nexttile; imshow(vizOverlay(cc(S.img), {cc(S.fluid1), cc(rm.mask & valid)}, ...
        {[0 1 0],[1 0 0]}, struct('fillAlpha',0.18,'lineWidth',2)));
    title(sprintf('modified (ours)   Dice %.2f', Em.dice));
end
figStyle(f,'fig_09_baseline_vs_modified');

%% ======================================================================
%  fig_10  Metrics comparison
%  ======================================================================
cmp = readtable(fullfile(p.tables,'comparison_vs_paper.csv'));
f = figure('Position',[20 20 1400 650]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Quantitative comparison on the Duke DME dataset (78 evaluable B-scans)');

nexttile;
M = [cmp.PaperReported, cmp.OurBaseline, cmp.OurModified];
b = bar(M); grid on
set(gca,'XTickLabel',cmp.Metric); ylabel('score'); ylim([0 1]);
legend({'Liu et al. 2021 (reported)','our baseline re-implementation','our modified pipeline'}, ...
    'Location','northoutside');
for k = 1:numel(b)
    text(b(k).XEndPoints, b(k).YEndPoints, compose('%.2f',b(k).YEndPoints), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
end
title('metrics');

nexttile;
tt = [mean(R.baseline.perScan.timeSec), mean(R.modified.perScan.timeSec), 1.2];
bar(tt,'FaceColor',[0.3 0.6 0.4]); grid on
set(gca,'XTickLabel',{'our baseline','our modified','paper (reported)'});
ylabel('seconds per B-scan');
text(1:3, tt, compose('%.2f s',tt), 'HorizontalAlignment','center','VerticalAlignment','bottom');
title('processing speed');
figStyle(f,'fig_10_metrics_comparison');

%% ======================================================================
%  fig_11  Per-scan Dice across the whole evaluable set
%  ======================================================================
ps = R.modified.perScan;
[~,ord] = sort(ps.diceA);
f = figure('Position',[20 20 1600 750]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Per-B-scan performance: one averaged number hides a lot of variance');

nexttile;
bar(ps.diceA(ord),'FaceColor',[0.25 0.45 0.75]); hold on
yline(mean(ps.diceA,'omitnan'),'r--','LineWidth',2);
grid on; xlabel('B-scan (sorted by Dice)'); ylabel('Dice vs expert A'); ylim([0 1]);
legend({'per-scan Dice',sprintf('mean %.3f',mean(ps.diceA,'omitnan'))},'Location','northwest');
title(sprintf('all %d evaluable B-scans', height(ps)));

nexttile;
scatter(ps.areaPx*0 + E.nFluid1, ps.diceA, 40, ps.subject, 'filled'); hold on
set(gca,'XScale','log'); grid on
xlabel('expert-A fluid area (pixels, log scale)'); ylabel('Dice');
c = colorbar; c.Label.String = 'subject';
ylim([0 1]);
title('performance scales with lesion size');
figStyle(f,'fig_11_per_scan_dice');

%% ======================================================================
%  fig_12  Inter-observer agreement -- the realistic ceiling
%  ======================================================================
% Two fellowship-trained graders do not agree with each other perfectly, so
% their mutual Dice is the practical upper bound for any automatic method
% scored against a single grader.  Neither paper reports this.
dAB = zeros(height(E),1);
for i = 1:height(E)
    S = loadBScan(E.subject(i), E.bscan(i));
    Ei = evalSegmentation(S.fluid2, S.fluid1, S.gradedMask);
    dAB(i) = Ei.dice;
end
f = figure('Position',[20 20 1300 620]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'How much do the two expert graders agree with each other?');

nexttile;
histogram(dAB,0:0.05:1,'FaceColor',[0.6 0.4 0.7]); hold on
xline(mean(dAB,'omitnan'),'k--','LineWidth',2);
xline(mean(ps.diceA,'omitnan'),'r--','LineWidth',2);
grid on; xlabel('Dice'); ylabel('number of B-scans');
legend({'expert A vs expert B', sprintf('their mean %.2f',mean(dAB,'omitnan')), ...
        sprintf('ours vs A %.2f',mean(ps.diceA,'omitnan'))},'Location','northwest');
title('inter-observer Dice');

nexttile;
scatter(dAB, ps.diceA, 45, 'filled','MarkerFaceAlpha',0.6); hold on
plot([0 1],[0 1],'k--'); grid on; axis square
xlabel('expert A vs expert B Dice'); ylabel('ours vs expert A Dice');
xlim([0 1]); ylim([0 1]);
title('we do better where the experts agree');
figStyle(f,'fig_12_interobserver');

save(fullfile(p.tables,'interobserver.mat'),'dAB');

%% ======================================================================
%  fig_13  Failure gallery -- shown honestly
%  ======================================================================
[~,worst] = sort(ps.diceA); worst = worst(1:4);
f = figure('Position',[10 10 1600 780]);
tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {'Failure cases, shown honestly (green = expert A, blue = expert B, red = ours)', ...
     'these are the four worst of the 78 evaluable B-scans'});
for i = 1:numel(worst)
    k = worst(i);
    S = loadBScan(ps.subject(k), ps.bscan(k));
    rm = segmentCME(S.img, cfgM, S.stack);
    valid = rm.roi.mask & S.gradedMask;
    a0 = max(1, round(min(rm.roi.ilm))-25);
    a1 = min(size(S.img,1), round(max(rm.roi.obm))+25);
    cc = @(X) X(a0:a1,:,:);
    nexttile;
    imshow(vizOverlay(cc(S.img), {cc(S.fluid1), cc(S.fluid2), cc(rm.mask & valid)}, ...
        {[0 1 0],[0.2 0.6 1],[1 0 0]}, struct('lineWidth',2)));
    title(sprintf('S%02d k%02d  Dice %.2f  (expert A marked %d px, expert B %d px)', ...
        ps.subject(k), ps.bscan(k), ps.diceA(k), nnz(S.fluid1), nnz(S.fluid2)));
end
figStyle(f,'fig_13_failure_gallery');

%% ======================================================================
%  fig_14  Why the CLAHE tile size matters
%  ======================================================================
S = loadBScan(5,24);
tiles = {[2 2],[4 4],[8 8],[16 16]};
f = figure('Position',[20 20 1700 860]);
tl = tiledlayout(f,2,4,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {'CLAHE tile size vs lesion size: with 8x8 tiles a tile is the same size as a cyst,', ...
     'so the local histogram equalisation flattens the very contrast we need'});
for i = 1:numel(tiles)
    c = cfgM; c.contrast.claheNumTiles = tiles{i};
    w = preprocessOCT(S.img, c);
    roi = extractROI(w, c);
    a0 = max(1,round(min(roi.ilm))-20); a1 = min(size(S.img,1),round(max(roi.obm))+20);
    nexttile(i);
    imshow(vizOverlay(w(a0:a1,:), {S.fluid1(a0:a1,:)}, {[0 1 0]}, struct('lineWidth',1)));
    tilePx = round([size(S.img,1) size(S.img,2)]./tiles{i});
    title(sprintf('%s tiles  (%dx%d px each)', mat2str(tiles{i}), tilePx(1), tilePx(2)));

    nexttile(i+4);
    histogram(w(S.fluid1 & roi.mask),0:6:255,'Normalization','probability'); hold on
    histogram(w(roi.mask & ~S.fluid1),0:6:255,'Normalization','probability');
    grid on; xlabel('gray'); ylim([0 0.12]);
    inV = w(S.fluid1 & roi.mask); outV = w(roi.mask & ~S.fluid1);
    sep = (mean(outV)-mean(inV))/sqrt(0.5*(var(inV)+var(outV)));
    title(sprintf('separation (CNR) = %.2f', sep));
    if i==1, legend({'cyst','retina'},'Location','northeast'); end
end
figStyle(f,'fig_14_preprocessing_effect');

fprintf('\nOur-contribution figures done.\n');
