%RUN_06_CLINICAL  The clinical output layer -- our addition, in neither paper.
%
%   fig_15_severity_distribution   how the 78 scans grade out
%   fig_16_severity_examples       example B-scans annotated with their numbers
%   fig_17_cyst_count_validation   our cyst count vs the experts' own count
%   table_04_clinical.csv          every per-scan clinical measurement
%
% The Duke manualFluid masks label each cyst with its own integer ID (values
% 1..15 rather than a plain 0/1 mask), so the graders effectively published a
% cyst COUNT as well as a cyst area.  That lets us validate the count output
% directly, which is a check neither paper could run.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

T = buildEvalSet(); E = T(T.evaluable,:);
L = load(fullfile(p.tables,'batchResults.mat')); R = L.R;
ps = R.modified.perScan;
cfg = cmeConfig('modified');

%% ----------------------------------------------------------------------
%  Recompute the data-driven grading cut-offs and state them explicitly
%  ----------------------------------------------------------------------
areas = ps.areaMm2(ps.areaMm2 > 0);
cuts  = prctile(areas, [33.3 66.7]);
fprintf(['Data-driven severity cut-offs (tertiles of the measured fluid area\n' ...
         'across the %d evaluable B-scans):  mild < %.4f <= moderate < %.4f <= severe  mm^2\n'], ...
         numel(areas), cuts(1), cuts(2));
fprintf('cmeConfig currently uses [%.2f %.2f] mm^2\n', cfg.clinical.areaCutMm2);

cfgT = cfg; cfgT.clinical.areaCutMm2 = cuts;

% regrade every scan with the recomputed cut-offs
grade = strings(height(ps),1);
for i = 1:height(ps)
    Ctmp.areaPx = ps.areaPx(i); Ctmp.areaMm2 = ps.areaMm2(i); Ctmp.cstUm = ps.cstUm(i);
    grade(i) = gradeSeverity(Ctmp, cfgT);
end
ps.gradeTertile = grade;
writetable(ps, fullfile(p.tables,'table_04_clinical.csv'));

%% ======================================================================
%  fig_15  Severity distribution
%  ======================================================================
f = figure('Position',[20 20 1500 650]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Automatic severity grading of the 78 evaluable B-scans (our clinical output layer)');

cats = ["none","mild","moderate","severe"];
cnt = arrayfun(@(g) sum(grade == g), cats);
nexttile;
b = bar(cnt,'FaceColor','flat');
b.CData = [0.7 0.7 0.7; 0.35 0.7 0.35; 0.95 0.75 0.2; 0.85 0.25 0.25];
set(gca,'XTickLabel',cats); ylabel('number of B-scans'); grid on
text(1:numel(cnt), cnt, string(cnt),'HorizontalAlignment','center','VerticalAlignment','bottom');
title('grade counts (area-based)');

nexttile;
histogram(areas, 25, 'FaceColor',[0.3 0.5 0.8]); hold on
xline(cuts(1),'g--','LineWidth',2); xline(cuts(2),'r--','LineWidth',2);
grid on; xlabel('segmented fluid area (mm^2)'); ylabel('B-scans');
legend({'measured area',sprintf('mild/moderate %.3f',cuts(1)), ...
        sprintf('moderate/severe %.3f',cuts(2))},'Location','northeast');
title('where the cut-offs sit');

nexttile;
scatter(ps.cstUm, ps.areaMm2, 45, 'filled','MarkerFaceAlpha',0.6); hold on
xline(250,'g--'); xline(400,'r--');
grid on; xlabel('central subfield thickness (\mum)'); ylabel('fluid area (mm^2)');
title(['area-based vs thickness-based grading' newline '(dashed = DRCR-style CST bands)']);
figStyle(f,'fig_15_severity_distribution');

%% ======================================================================
%  fig_16  Annotated examples, one per grade
%  ======================================================================
f = figure('Position',[10 10 1600 900]);
tl = tiledlayout(f,3,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Example B-scans with the automatic clinical read-out printed on the image');

shown = 0;
tileIdx = [1 3 5];          % left column of the 3x2 layout
for g = ["mild","moderate","severe"]
    idx = find(grade == g);
    if isempty(idx), continue; end
    [~,mid] = min(abs(ps.areaMm2(idx) - median(ps.areaMm2(idx))));
    k = idx(mid);
    S = loadBScan(ps.subject(k), ps.bscan(k));
    rm = segmentCME(S.img, cfg, S.stack);
    valid = rm.roi.mask & S.gradedMask;
    a0 = max(1,round(min(rm.roi.ilm))-25); a1 = min(size(S.img,1),round(max(rm.roi.obm))+25);
    cc = @(X) X(a0:a1,:,:);
    shown = shown + 1;
    nexttile(tileIdx(shown));
    imshow(vizOverlay(cc(S.img), {cc(S.fluid1), cc(rm.mask & valid)}, {[0 1 0],[1 0 0]}, ...
        struct('fillAlpha',0.22,'lineWidth',2)));
    txt = sprintf('GRADE: %s\narea %.3f mm^2\ncysts %d\nmean diam %.0f um\nCST %.0f um', ...
        upper(g), ps.areaMm2(k), ps.cystCount(k), ps.meanDiamUm(k), ps.cstUm(k));
    text(10, 18, txt, 'Color','y','FontSize',11,'FontName','Consolas', ...
        'VerticalAlignment','top','BackgroundColor',[0 0 0 ],'Margin',3);
    title(sprintf('S%02d k%02d  --  %s', ps.subject(k), ps.bscan(k), upper(g)));
end

% the grading rule itself, written out, spanning the whole right column
nexttile(2,[3 1]); axis off
txt = sprintf([ ...
 'HOW THE GRADE IS COMPUTED\n\n' ...
 'AREA-BASED (project-defined, data-driven)\n' ...
 '  mild      area <  %.3f mm^2\n' ...
 '  moderate  %.3f <= area < %.3f mm^2\n' ...
 '  severe    area >= %.3f mm^2\n' ...
 '  Cut-offs are the TERTILES of the fluid area we\n' ...
 '  measured across the 78 evaluable Duke B-scans.\n' ...
 '  There is no consensus clinical cut-off for cyst\n' ...
 '  AREA on a single B-scan, so we do not claim one.\n\n' ...
 'CST-BASED (literature anchored)\n' ...
 '  mild      CST <  250 um\n' ...
 '  moderate  250 <= CST < 400 um\n' ...
 '  severe    CST >= 400 um\n' ...
 '  Central subfield thickness is what commercial OCT\n' ...
 '  software reports and what DRCR.net Protocol T uses\n' ...
 '  to define centre-involved DME.\n\n' ...
 'Area is converted with the Duke pixel pitch:\n' ...
 '  %.2f um axial x %.2f um lateral = %.2e mm^2/px\n' ...
 'Equivalent-circle diameter is derived from the\n' ...
 'PHYSICAL area, so the 3:1 pixel anisotropy is handled.'], ...
 cuts(1), cuts(1), cuts(2), cuts(2), ...
 cfg.pixelPitchAxialUm, cfg.pixelPitchLateralUm, cfg.pixelAreaMm2);
text(0.02,0.98,txt,'Units','normalized','VerticalAlignment','top', ...
    'FontName','Consolas','FontSize',10);
figStyle(f,'fig_16_severity_examples');

%% ======================================================================
%  fig_17  Cyst count validation against the graders' own instance labels
%  ======================================================================
gtCount = zeros(height(E),1);
for i = 1:height(E)
    m = matfile(fullfile(p.data, sprintf('Subject_%02d.mat', E.subject(i))));
    F = m.manualFluid1(:,:,E.bscan(i));
    gtCount(i) = numel(unique(F(F > 0)));     % the graders numbered each cyst
end

f = figure('Position',[20 20 1300 620]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Cyst count: our connected-component count vs the graders'' own cyst IDs');

nexttile;
scatter(gtCount, ps.cystCount, 55,'filled','MarkerFaceAlpha',0.6); hold on
mx = max([gtCount; ps.cystCount])+1;
plot([0 mx],[0 mx],'k--'); grid on; axis square
xlabel('expert A cyst count'); ylabel('our cyst count');
cc = corr(gtCount, ps.cystCount, 'rows','complete');
% Report what this actually shows.  Connected-component counting is not the
% same operation the graders performed: they split a confluent cyst cluster by
% eye into separate cysts, whereas a connected-component labeller merges
% anything that touches.  The count output is therefore NOT reliable, and we
% say so instead of leaving a scatter that looks like agreement.
title(sprintf(['CYST COUNT DOES NOT AGREE  (Pearson r = %.2f)' newline ...
    'component labelling splits confluent clusters differently from a human'], cc));

nexttile;
scatter(ps.areaMm2, arrayfun(@(i) nnz(loadFluidArea(p, E.subject(i), E.bscan(i)))*cfg.pixelAreaMm2, ...
    (1:height(E))'), 55, 'filled','MarkerFaceAlpha',0.6); hold on
mx = max(ps.areaMm2)*1.1;
plot([0 mx],[0 mx],'k--'); grid on; axis square
xlabel('our fluid area (mm^2)'); ylabel('expert A fluid area (mm^2)');
gtArea = arrayfun(@(i) nnz(loadFluidArea(p, E.subject(i), E.bscan(i)))*cfg.pixelAreaMm2, (1:height(E))');
ca = corr(ps.areaMm2, gtArea, 'rows','complete');
title(sprintf(['AREA does agree  (Pearson r = %.2f)' newline ...
    'points below the line = we over-segment relative to expert A'], ca));
figStyle(f,'fig_17_cyst_count_validation');

writetable(table(E.subject, E.bscan, gtCount, ps.cystCount, ...
    'VariableNames',{'subject','bscan','expertCystCount','ourCystCount'}), ...
    fullfile(p.tables,'table_05_cyst_counts.csv'));

fprintf('\nClinical figures done.\n');


function F = loadFluidArea(p, s, k)
m = matfile(fullfile(p.data, sprintf('Subject_%02d.mat', s)));
F = m.manualFluid1(:,:,k) > 0;
end
