%RUN_03_FIGURES_PAPER  Generate the figures that mirror the two reference
%                      papers, using OUR data and OUR results.
%
% Produces:
%   fig_01_pipeline_overview      stage-by-stage walk through the pipeline
%   fig_02_wave_maps              seawater / wave / coastline maps (2020 Fig.2)
%   fig_03_roi_extraction         row-mean-gray curve + A/B cuts (2021 Fig.3)
%   fig_04_directional_contours   the 4 directions + fusion (2021 Fig.5)
%   fig_05_segmentation_gallery   final overlays vs both experts, many subjects
%   fig_06_roi_vs_expert_layers   our ILM/OBM against the expert traces

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

DEMO = [5 24];                      % the representative scan used throughout
S = loadBScan(DEMO(1), DEMO(2));
res = segmentCME(S.img, cfg, S.stack);

% Crop rows to the retina so every picture is readable in a slide.
r0 = max(1, round(min(res.roi.ilm)) - 30);
r1 = min(size(S.img,1), round(max(res.roi.obm)) + 30);
cr = @(X) X(r0:r1, :, :);

%% ======================================================================
%  fig_01  Pipeline overview, stage by stage
%  ======================================================================
f = figure('Position',[20 20 1700 760]);
tl = tiledlayout(f,3,3,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, sprintf('CME segmentation pipeline, stage by stage (Subject %02d, B-scan %02d)', DEMO));

nexttile; imshow(cr(res.stages.raw),[0 255]);      title('1. raw B-scan');
nexttile; imshow(cr(res.stages.denoised),[0 255]); title('2. denoised (median 5x5)');
nexttile; imshow(cr(res.stages.enhanced),[0 255]); title('3. contrast enhanced (CLAHE)');

nexttile;
imshow(vizOverlay(cr(res.work), {cr(res.roi.mask), cr(res.roi.innerBand)}, ...
    {[0 0.8 1],[1 0.9 0]}, struct('fillAlpha',0.15)));
title('4. ROI: retina (cyan) and inner band (yellow)');

nexttile;
imshow(vizOverlay(cr(res.work), {cr(imdilate(res.dirs.fused,strel('disk',1)))}, {[1 0.85 0.1]}, ...
    struct('fillAlpha',1)));
title('5. fused 4-direction wave contours');

nexttile;
imshow(vizOverlay(cr(res.work), {cr(res.integrate.candidates)}, {[0 1 1]}, ...
    struct('fillAlpha',0.45,'lineWidth',1)));
title('6. candidate basins (before screening)');

nexttile;
imshow(vizOverlay(cr(res.work), {cr(res.mask)}, {[1 0.3 0.3]}, ...
    struct('fillAlpha',0.45,'lineWidth',1)));
title('7. after gray / size / contrast screening');

nexttile;
imshow(vizOverlay(cr(S.img), {cr(S.fluid1), cr(res.mask)}, {[0 1 0],[1 0 0]}, ...
    struct('fillAlpha',0.20,'lineWidth',2)));
title('8. final overlay: expert A (green), ours (red)');

nexttile; axis off;
C = res.clinical;
txt = sprintf(['CLINICAL OUTPUT\n\n' ...
    'fluid area   %.0f px  (%.3f mm^2)\n' ...
    'cyst count   %d\n' ...
    'mean diam.   %.0f um\n' ...
    'max diam.    %.0f um\n' ...
    'CST          %.0f um\n' ...
    'retina thick %.0f um\n\n' ...
    'GRADE:  %s\n(%s)'], ...
    C.areaPx, C.areaMm2, C.cystCount, C.meanDiamUm, C.maxDiamUm, ...
    C.cstUm, C.retinalThicknessUm, upper(C.grade), C.gradeBasis);
text(0.02,0.98,txt,'Units','normalized','VerticalAlignment','top', ...
    'FontName','Consolas','FontSize',11);
figStyle(f,'fig_01_pipeline_overview');

%% ======================================================================
%  fig_02  Wave maps: seawater / wave / coastline  (mirrors Lou 2020 Fig. 2)
%  ======================================================================
W = res.dirs.exemplar;              % the theta = 0 sweep
f = figure('Position',[20 20 1700 660]);
tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Wave potential energy and the correction equation (\theta = 0, sweeping downwards)');

nexttile; imshow(cr(res.work),[0 255]); title('(a) working image');

nexttile;
imagesc(cr(W.wave)); axis image off; colormap(gca,parula); colorbar;
title('(b) wave potential energy  \phi_v + \phi_g');

% The three-region "wave map" of the 2020 paper, colour-coded the same way:
% blue = seawater (noise), green = wave (boundary area), red = coastline.
wm = zeros([size(cr(W.IT)) 3]);
sea   = cr(~W.IT & ~W.Pen);
wave  = cr(W.IT);
coast = cr(W.Pen & ~W.IT);
wm(:,:,3) = double(sea);
wm(:,:,2) = double(wave);
wm(:,:,1) = double(coast);
nexttile; imshow(wm);
title('(c) wave map: seawater(b) / wave(g) / coastline(r)');

nexttile; imagesc(cr(W.phi_v)); axis image off; colorbar;
title('(d) kinetic energy  \phi_v = v\cdotv_q\cdot\sigma');
nexttile; imagesc(cr(W.sigma)); axis image off; colorbar;
title('(e) regulation factor \sigma');

nexttile;
imshow(vizOverlay(cr(res.work), {cr(imdilate(W.bnd,strel('disk',1)))}, {[1 0.2 0.2]}, struct('fillAlpha',1)));
title('(f) boundary after the correction equation');
figStyle(f,'fig_02_wave_maps');

%% ======================================================================
%  fig_03  ROI extraction  (mirrors Liu 2021 Fig. 3)
%  ======================================================================
f = figure('Position',[20 20 1500 620]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'ROI extraction from the row-mean-gray distribution (Eq. 1)');

nexttile;
plot(res.roi.rowMeanRaw, 1:numel(res.roi.rowMeanRaw), 'Color',[.7 .7 .7]); hold on
plot(res.roi.rowMean,    1:numel(res.roi.rowMean), 'k','LineWidth',1.5);
xline(res.roi.imgMean,'b','LineWidth',2);
yline(res.roi.A,'r--','LineWidth',1.5); yline(res.roi.B,'r--','LineWidth',1.5);
text(res.roi.imgMean, res.roi.A-8, '  A', 'Color','r','FontWeight','bold');
text(res.roi.imgMean, res.roi.B+14,'  B', 'Color','r','FontWeight','bold');
set(gca,'YDir','reverse'); ylim([1 numel(res.roi.rowMean)]); grid on
xlabel('mean gray of the row'); ylabel('row (depth)');
legend({'raw row mean','smoothed','whole-image mean','A / B'},'Location','southeast');
title('row-mean-gray curve');

nexttile;
imshow(S.img,[0 255]); hold on
yline(res.roi.A,'r--','LineWidth',1.5); yline(res.roi.B,'r--','LineWidth',1.5);
title('A and B on the B-scan');

nexttile;
imshow(vizOverlay(S.img, {res.roi.mask}, {[0 0.9 1]}, struct('fillAlpha',0.18))); hold on
% Keep the line handles and pass them to legend: imshow puts an image object
% into the axes first, and legend would otherwise label that instead.
hI = plot(res.roi.ilm,'g','LineWidth',1.5);
hO = plot(res.roi.obm,'y','LineWidth',1.5);
legend([hI hO], {'ILM','OBM'}, 'TextColor','w','Color',[0 0 0],'Location','southwest');
title('refined ROI (ILM to OBM)');
figStyle(f,'fig_03_roi_extraction');

%% ======================================================================
%  fig_04  Directional contours + fusion  (mirrors Liu 2021 Fig. 5)
%  ======================================================================
f = figure('Position',[20 20 1700 660]);
tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Omnidirectional wave operator: four operating directions and their fusion');
lbl = {'(a) \theta = 0','(b) \theta = \pi/2','(c) \theta = \pi','(d) \theta = 3\pi/2'};
for t = 1:4
    nexttile;
    imshow(vizOverlay(cr(res.work), {cr(res.dirs.dirBnd{t})}, {[1 0.25 0.25]}, struct('lineWidth',1)));
    title(sprintf('%s  (%.1f%% of ROI)', lbl{t}, 100*mean(res.dirs.dirBnd{t}(res.roi.mask))));
end
nexttile;
imshow(vizOverlay(cr(res.work), {cr(res.dirs.fused)}, {[1 0.9 0.1]}, struct('lineWidth',1)));
title(sprintf('(e) fusion of all four  (%.1f%%)', 100*mean(res.dirs.fused(res.roi.mask))));
nexttile;
imshow(vizOverlay(cr(S.img), {cr(S.fluid1), cr(res.mask)}, {[0 1 0],[1 0 0]}, ...
    struct('fillAlpha',0.20,'lineWidth',2)));
title('(f) final result vs expert A');
figStyle(f,'fig_04_directional_contours');

%% ======================================================================
%  fig_05  Segmentation gallery across subjects and severities
%  ======================================================================
T = buildEvalSet(); E = T(T.evaluable,:);
% pick 8 scans spanning the range of fluid burden, from different subjects
[~,ord] = sort(E.nFluid1);
pick = E(ord(round(linspace(1, height(E), 8))), :);

f = figure('Position',[10 10 1700 1150]);
tl = tiledlayout(f,4,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Final segmentation vs expert graders (green = expert A, blue = expert B, red = ours)');
for i = 1:height(pick)
    Si = loadBScan(pick.subject(i), pick.bscan(i));
    ri = segmentCME(Si.img, cfg, Si.stack);
    valid = ri.roi.mask & Si.gradedMask;
    Ei = evalSegmentation(ri.mask & valid, Si.fluid1, valid);
    a0 = max(1, round(min(ri.roi.ilm))-25);
    a1 = min(size(Si.img,1), round(max(ri.roi.obm))+25);
    cc = @(X) X(a0:a1,:,:);
    nexttile;
    imshow(vizOverlay(cc(Si.img), {cc(Si.fluid1), cc(Si.fluid2), cc(ri.mask & valid)}, ...
        {[0 1 0],[0.2 0.6 1],[1 0 0]}, struct('lineWidth',2)));
    title(sprintf('S%02d k%02d   Dice %.2f  P %.2f  R %.2f', ...
        pick.subject(i), pick.bscan(i), Ei.dice, Ei.precision, Ei.recall));
end
figStyle(f,'fig_05_segmentation_gallery');

%% ======================================================================
%  fig_06  ROI accuracy against the expert layer traces
%  ======================================================================
LT = T(logical(T.hasLayers),:);
idx = unique(max(1, min(height(LT), round(linspace(1,height(LT),4)))));
sel = LT(idx,:);
f = figure('Position',[20 20 1500 720]);
tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'ROI boundaries vs expert layer traces (solid = ours, dashed = expert)');
for i = 1:height(sel)
    Si = loadBScan(sel.subject(i), sel.bscan(i));
    wi = preprocessOCT(Si.img, cfg); ri = extractROI(wi, cfg);
    gI = Si.layers1(1,:); gB = Si.layers1(8,:);
    ok = isfinite(gI) & isfinite(gB);
    nexttile; imshow(Si.img,[0 255]); hold on
    plot(ri.ilm,'g','LineWidth',1.5); plot(ri.obm,'y','LineWidth',1.5);
    plot(gI,'g--','LineWidth',1.2);   plot(gB,'y--','LineWidth',1.2);
    title(sprintf('S%02d k%02d  ILM err %.1f px, OBM err %.1f px', ...
        sel.subject(i), sel.bscan(i), ...
        mean(abs(ri.ilm(ok)-gI(ok))), mean(abs(ri.obm(ok)-gB(ok)))));
end
figStyle(f,'fig_06_roi_vs_expert_layers');

fprintf('\nPaper-mimicking figures done.\n');
