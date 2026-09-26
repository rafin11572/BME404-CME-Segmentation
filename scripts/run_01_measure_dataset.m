%RUN_01_MEASURE_DATASET  Measure the properties of the expert annotations that
%                        the pipeline's screening rules depend on.
%
% These are the numbers quoted in docs/FINDINGS.md sections 3 and 4, and the
% evidence behind our two departures from the published parameter values.
% Saves cystThicknessStats.mat, which run_04_ablation needs.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

T = buildEvalSet(); E = T(T.evaluable,:);

rel = []; areas = []; depth = []; heights = [];
fprintf('measuring %d evaluable B-scans ...\n', height(E));
for i = 1:height(E)
    S = loadBScan(E.subject(i), E.bscan(i));
    work = preprocessOCT(S.img, cfg);
    roi  = extractROI(work, cfg);

    st = regionprops(S.fluid1, 'BoundingBox','Area','Centroid');
    for k = 1:numel(st)
        rel(end+1)     = st(k).BoundingBox(4) / roi.thicknessPx;  %#ok<SAGROW>
        areas(end+1)   = st(k).Area;                              %#ok<SAGROW>
        heights(end+1) = st(k).BoundingBox(4);                    %#ok<SAGROW>
    end

    % relative depth of every annotated fluid PIXEL
    [rr,cc] = find(S.fluid1 & roi.mask);
    if ~isempty(rr)
        d = (rr - roi.ilm(cc)') ./ max(roi.obm(cc)' - roi.ilm(cc)', 1);
        depth = [depth; d]; %#ok<AGROW>
    end
    if mod(i,20)==0, fprintf('  %d/%d\n', i, height(E)); end
end

fprintf('\n--- cyst size, relative to mean retinal thickness (%d cysts) ---\n', numel(rel));
fprintf('  p5 %.3f  p25 %.3f  median %.3f  p75 %.3f  p95 %.3f  max %.3f\n', ...
    prctile(rel,[5 25 50 75 95 100]));
fprintf('  inside the paper''s [0.5 1.5] window : %.1f%%\n', 100*mean(rel>=0.5 & rel<=1.5));
fprintf('  inside our        [0.08 1.2] window : %.1f%%\n', 100*mean(rel>=0.08 & rel<=1.2));

fprintf('\n--- fluid depth inside the retina (0 = ILM, 1 = OBM) ---\n');
fprintf('  p5 %.2f  p25 %.2f  median %.2f  p75 %.2f  p95 %.2f\n', prctile(depth,[5 25 50 75 95]));
fprintf('  our inner band is [%.2f %.2f] -> captures %.1f%% of fluid pixels\n', ...
    cfg.roi.depthBand, 100*mean(depth>=cfg.roi.depthBand(1) & depth<=cfg.roi.depthBand(2)));

fprintf('\n--- cyst area ---\n');
fprintf('  median %.0f px (%.4f mm^2), p95 %.0f px\n', ...
    median(areas), median(areas)*cfg.pixelAreaMm2, prctile(areas,95));

save(fullfile(p.tables,'cystThicknessStats.mat'), 'rel','areas','depth','heights');
writetable(table(rel(:), areas(:), heights(:), ...
    'VariableNames',{'thicknessRatio','areaPx','heightPx'}), ...
    fullfile(p.tables,'table_01_expert_cyst_stats.csv'));
fprintf('\nSaved cystThicknessStats.mat\n');
