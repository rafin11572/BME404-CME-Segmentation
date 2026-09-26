%DEV_CEILING  What is the best Dice ANY global intensity threshold could give?
% DEVELOPMENT / MEASUREMENT script -- not part of the deliverable pipeline.
% It produced numbers quoted in docs/FINDINGS.md; kept so they are reproducible.
% If the oracle threshold is already poor, the problem is not the threshold --
% it is that cysts and other dark tissue overlap in intensity, and the method
% needs spatial context, not a better cut point.
clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

devSet = [1 31; 4 22; 5 24; 6 26; 8 24; 10 21; 9 26; 2 29];
base = cmeConfig('modified');
base.contrast.claheNumTiles = [4 4];

fprintf('%-10s %8s %8s %8s %9s %8s %8s\n', ...
    'scan','ROIpx','GTpx','GT%ROI','bestDice','bestPct','AUC');
allIn = []; allOut = [];
for i = 1:size(devSet,1)
    S = loadBScan(devSet(i,1), devSet(i,2));
    work = preprocessOCT(S.img, base);
    roi  = extractROI(work, base);
    gt   = S.fluid1 & roi.mask;

    g = work(roi.mask);
    inV  = work(gt);
    outV = work(roi.mask & ~gt);
    allIn = [allIn; inV]; allOut = [allOut; outV]; %#ok<AGROW>

    best = 0; bestP = 0;
    for pct = 1:60
        lvl = prctile(g, pct);
        M = (work < lvl) & roi.mask;
        E = evalSegmentation(M, gt, roi.mask);
        if E.dice > best, best = E.dice; bestP = pct; end
    end
    % separability: P(random cyst pixel darker than random non-cyst pixel)
    auc = mean(randsample(inV,min(3000,numel(inV))) < ...
               randsample(outV,min(3000,numel(inV)))');
    fprintf('S%02d k%02d %8d %8d %7.1f%% %9.3f %8d %8.3f\n', ...
        devSet(i,1), devSet(i,2), nnz(roi.mask), nnz(gt), ...
        100*nnz(gt)/nnz(roi.mask), best, bestP, mean(auc(:)));
end

fprintf('\npooled: cyst gray mean %.1f (sd %.1f), non-cyst %.1f (sd %.1f)\n', ...
    mean(allIn), std(allIn), mean(allOut), std(allOut));

% picture of what competes with the cysts
S = loadBScan(5,24);
work = preprocessOCT(S.img, base); roi = extractROI(work, base);
g = work(roi.mask); lvl = prctile(g,15);
f = figure('Color','w','Position',[20 20 1600 900]);
tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
nexttile; imshow(work,[0 255]); hold on; visboundaries(S.fluid1,'Color','g','LineWidth',0.8); title('work + expert A');
nexttile; imshow((work<lvl) & roi.mask); title('darkest 15% inside ROI');
nexttile; imshow(roi.mask); title('ROI mask');
nexttile;
histogram(work(S.fluid1 & roi.mask),0:4:255,'Normalization','probability'); hold on
histogram(work(roi.mask & ~S.fluid1),0:4:255,'Normalization','probability');
legend({'cyst','retina non-cyst'}); xlabel('gray'); title('intensity overlap'); grid on
exportgraphics(f, fullfile(p.figures,'dev_ceiling.png'),'Resolution',140);
disp('done');
