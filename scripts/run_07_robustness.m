%RUN_07_ROBUSTNESS  How often does the pipeline invent fluid on a healthy scan?
%
% Neither reference paper reports this, but it is the number a screening tool
% would actually be judged on: on a B-scan with no expert-marked fluid, the
% pipeline should output nothing.  We run it on a random sample of the 532
% fluid-free B-scans and report the false-positive rate.
%
% Outputs: fig_18_false_positive_rate.png, table_06_robustness.csv

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

T = buildEvalSet();
% CRITICAL: a B-scan with no fluid in manualFluid1 is NOT necessarily
% fluid-free.  Only 110 of the 610 Duke B-scans were graded at all; on the
% other 500 the fluid mask is entirely NaN.  The genuine negatives are the
% B-scans the experts DID examine and found no fluid in, identified by having
% finite layer traces but no fluid from either
% grader.  There are 19 such B-scans once BOTH graders are taken into account
% (110 graded, 78 with expert-A fluid, 86 with expert-B fluid).  Scoring
% against ungraded scans would be meaningless.
neg = T(logical(T.hasLayers) & ~T.hasFluid1 & ~T.hasFluid2, :);
nSample = height(neg);
sel = neg;

fprintf('testing %d fluid-free B-scans ...\n', nSample);
fpPx = zeros(nSample,1); fpFrac = zeros(nSample,1); nComp = zeros(nSample,1);
for i = 1:nSample
    S = loadBScan(sel.subject(i), sel.bscan(i));
    r = segmentCME(S.img, cfg, S.stack);
    m = r.mask & r.roi.mask & S.gradedMask;
    fpPx(i)   = nnz(m);
    fpFrac(i) = nnz(m) / max(nnz(r.roi.mask & S.gradedMask),1);
    nComp(i)  = r.clinical.cystCount;
    if mod(i,20)==0, fprintf('  %d/%d\n', i, nSample); end
end

cleanRate = mean(fpPx == 0);
smallRate = mean(fpPx < 200);
fprintf('\nFluid-free B-scans:\n');
fprintf('  completely clean (0 px)        : %.1f%%\n', 100*cleanRate);
fprintf('  under 200 px of false fluid    : %.1f%%\n', 100*smallRate);
fprintf('  median false-positive area     : %.0f px (%.4f mm^2)\n', ...
    median(fpPx), median(fpPx)*cfg.pixelAreaMm2);
fprintf('  mean false fluid as %% of retina: %.2f%%\n', 100*mean(fpFrac));

Tab = table(sel.subject, sel.bscan, fpPx, fpFrac, nComp, ...
    'VariableNames',{'subject','bscan','falsePositivePx','falseFractionOfROI','components'});
writetable(Tab, fullfile(p.tables,'table_06_robustness.csv'));

f = figure('Position',[20 20 1400 620]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {sprintf('Robustness on the %d B-scans the experts GRADED and found fluid-free', nSample), ...
     'the other 500 Duke B-scans were never graded, so they cannot be scored at all'});

nexttile;
histogram(fpPx, 30, 'FaceColor',[0.35 0.6 0.4]); hold on
xline(median(fpPx),'r--','LineWidth',2);
grid on; xlabel('false-positive area (pixels)'); ylabel('B-scans');
legend({'false fluid area', sprintf('median %.0f px', median(fpPx))});
title(sprintf('%.0f%% produce under 200 px of false fluid', 100*smallRate));

nexttile;
% Compare against the true-positive areas so the scale is meaningful.
L = load(fullfile(p.tables,'batchResults.mat')); ps = L.R.modified.perScan;
histogram(ps.areaPx, 30, 'Normalization','probability','FaceColor',[0.8 0.3 0.3]); hold on
histogram(fpPx,      30, 'Normalization','probability','FaceColor',[0.35 0.6 0.4]);
set(gca,'XScale','log'); grid on
xlabel('segmented area (px, log scale)'); ylabel('fraction of B-scans');
legend({'scans WITH expert fluid','scans WITHOUT expert fluid'},'Location','northwest');
% Say what the data actually shows, not what we would like it to show.
ovl = 1 - 0.5*sum(abs(histcounts(log10(max(ps.areaPx,1)),20,'Normalization','probability') - ...
                      histcounts(log10(max(fpPx,1)),   20,'Normalization','probability')));
title(sprintf(['the two populations OVERLAP heavily (%.0f%%)' newline ...
    'segmented area alone cannot tell a healthy B-scan from a diseased one'], 100*ovl));
figStyle(f,'fig_18_false_positive_rate');

fprintf('\nRobustness check done.\n');
