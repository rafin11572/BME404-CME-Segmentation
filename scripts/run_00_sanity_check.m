%RUN_00_SANITY_CHECK  Confirm the data loads, the pairing is right, and the
%                     pipeline runs end-to-end on one B-scan.
%
% Run this first.  It answers three questions:
%   1. Do the .mat files load and does images(:,:,k) really correspond to
%      manualFluid1(:,:,k)?  (checked visually -- the mask must sit on dark
%      cystoid spaces inside the retina, not floating in the vitreous)
%   2. Does the ROI extraction find the retina?
%   3. Does the wave operator produce plausible contours?

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

subject = 5;  bscan = 24;     % a B-scan known to carry expert-A fluid
S = loadBScan(subject, bscan);

fprintf('Subject %02d, B-scan %02d\n', subject, bscan);
fprintf('  image   : %s, range [%g %g]\n', mat2str(size(S.img)), min(S.img(:)), max(S.img(:)));
fprintf('  fluid A : %d px    fluid B : %d px\n', nnz(S.fluid1), nnz(S.fluid2));
fprintf('  layers  : %s\n', string(S.hasLayers));

cfg = cmeConfig('modified');
res = segmentCME(S.img, cfg, S.stack);

fprintf('\nPipeline timing (s):\n');
disp(res.timing);
fprintf('ROI: A=%d B=%d, mean retinal thickness = %.1f px (%.0f um)\n', ...
    res.roi.A, res.roi.B, res.roi.thicknessPx, res.roi.thicknessPx*cfg.pixelPitchAxialUm);
fprintf('Segmented %d px in %d components\n', nnz(res.mask), res.clinical.cystCount);

% Score inside the retina AND inside the window the graders actually
% annotated -- see loadBScan for why the second restriction is necessary.
valid = res.roi.mask & S.gradedMask;
E  = evalSegmentation(res.mask & valid, S.fluid1, valid);
EB = evalSegmentation(res.mask & valid, S.fluid2, valid);
EX = evalSegmentation(S.fluid2, S.fluid1, valid);
fprintf('\nvs Expert A : Dice %.3f  Prec %.3f  Rec %.3f  Acc %.3f\n', ...
    E.dice, E.precision, E.recall, E.accuracy);
fprintf('vs Expert B : Dice %.3f  Prec %.3f  Rec %.3f\n', EB.dice, EB.precision, EB.recall);
fprintf('A vs B (inter-observer ceiling): Dice %.3f\n', EX.dice);

% ---- quick visual -----------------------------------------------------
f = figure('Color','w','Position',[80 80 1500 820]);
tiledlayout(f,2,3,'TileSpacing','compact','Padding','loose');

nexttile; imshow(S.img,[0 255]); title('raw B-scan');
nexttile; imshow(res.work,[0 255]); title('working image (median + CLAHE)');

nexttile; imshow(S.img,[0 255]); hold on;
plot(res.roi.ilm,'g','LineWidth',1.5); plot(res.roi.obm,'c','LineWidth',1.5);
yline(res.roi.A,'y--'); yline(res.roi.B,'y--');
title('ROI: ILM (green) / OBM (cyan)');

nexttile; imshow(res.dirs.fused); title('fused 4-direction contours');
nexttile; imshow(res.mask);       title('final mask after screening');

nexttile; imshow(S.img,[0 255]); hold on;
visboundaries(S.fluid1,'Color','g','LineWidth',1);
visboundaries(res.mask, 'Color','r','LineWidth',1);
title('green = expert A, red = ours');

exportgraphics(f, fullfile(p.figures,'fig_00_sanity_check.png'), 'Resolution',200);
fprintf('\nSaved %s\n', fullfile(p.figures,'fig_00_sanity_check.png'));
