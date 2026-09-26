%RUN_P_PRESENTATION  Build the slide-deck figure set.
%
%   Writes results/figures/presentation/P01..P08.png
%
% WHY A SEPARATE SET
%   The figures in results/figures/ are the technical record: wave maps,
%   energy terms, ablations, protocol analyses.  They belong in the report,
%   where there is room to explain them.  A live demo to an audience that has
%   never seen an OCT B-scan needs something different -- pictures that carry
%   their own meaning in five seconds, and at most one number per slide.
%
% RULES THIS SET FOLLOWS
%   1. Our modified pipeline only.  No baseline re-implementation curves.
%   2. Pictures before plots.  The most convincing evidence for a lay audience
%      is "here is the scan, here is what the doctor drew, here is what our
%      code found, look how they overlap".
%   3. At most one idea per figure, and no equation, colorbar, log axis or
%      ROC-style plot anywhere.
%   4. Every number shown is the honest one from the full evaluation -- the
%      framing is simplified, the measurements are not.
%
% Everything here reads the same results/tables/ files as the technical set,
% so the two can never disagree.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
outDir = fullfile(p.figures, 'presentation');
if ~exist(outDir,'dir'), mkdir(outDir); end

cfg = cmeConfig('modified');
T   = buildEvalSet();
E   = T(T.evaluable,:);
L   = load(fullfile(p.tables,'batchResults.mat'));  ps = L.R.modified.perScan;

% Colours used consistently across the whole deck.
COL.expert = [0.15 0.85 0.15];
COL.ours   = [1.00 0.20 0.20];
COL.bar    = [0.20 0.45 0.75];
COL.accent = [0.90 0.55 0.15];

save_ = @(f,name) exportgraphics(f, fullfile(outDir,[name '.png']), ...
    'Resolution',200, 'BackgroundColor','w');


%% =====================================================================
%  P01  What are we looking at, and what is the disease?
%  =====================================================================
% Opens the talk.  Nobody in the room has necessarily seen an OCT B-scan.
Sd = loadBScan(6, 33, 1);                       % clear, obvious cysts
Sn = loadBScan(6, 5, 1);                        % same eye, a quiet slice
rd = segmentCME(Sd.img, cfg, Sd.stack);
[a0,a1,ac0,ac1] = tightCrop(rd, Sd.fluid1);

f = figure('Color','w','Position',[40 40 1500 530]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'A retinal OCT scan, and what cystoid macular edema looks like');

nexttile; imshow(Sn.img(a0:a1, ac0:ac1),[0 255]);
title('healthy-looking retina: smooth, evenly layered','FontSize',12);

nexttile; imshow(Sd.img(a0:a1, ac0:ac1),[0 255]); hold on
st = regionprops(Sd.fluid1(a0:a1, ac0:ac1), 'Centroid','Area');
[~,big] = max([st.Area]);
cxy = st(big).Centroid;
annotation_arrow(cxy, 'fluid pockets (cysts)');
title('the same eye with edema: dark fluid pockets','FontSize',12);
save_(f,'P01_what_is_cme');
fprintf('  P01 done\n');


%% =====================================================================
%  P02  How the method works, told in five pictures and no equations
%  =====================================================================
S = loadBScan(5, 24, 1);
r = segmentCME(S.img, cfg, S.stack);
% Crop BOTH ways.  A full B-scan is about 4:1, so in a five-across layout each
% panel collapses to a thin strip with most of the canvas empty.  Cropping
% laterally around the lesion makes each panel roughly 2:1 and fills the slide.
[b0,b1,bc0,bc1] = tightCrop(r, S.fluid1);
cr = @(X) X(b0:b1, bc0:bc1, :);

f = figure('Color','w','Position',[20 20 1700 330]);
tl = tiledlayout(f,1,5,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'How our method finds the fluid, step by step');

nexttile; imshow(cr(r.stages.raw),[0 255]);
title({'1. the scan','as the machine gives it'},'FontSize',11);

nexttile; imshow(cr(r.stages.denoised),[0 255]);
title({'2. clean up the speckle','(median filter)'},'FontSize',11);

nexttile;
imshow(vizOverlay(cr(r.work), {cr(r.roi.innerBand)}, {[0.2 0.7 1]}, ...
    struct('fillAlpha',0.30)));
title({'3. find the retina','ignore everything else'},'FontSize',11);

nexttile;
imshow(vizOverlay(cr(r.work), {cr(r.mask)}, {COL.ours}, ...
    struct('fillAlpha',0.55,'lineWidth',1)));
title({'4. pick out the dark','fluid pockets'},'FontSize',11);

nexttile;
imshow(vizOverlay(cr(S.img), {cr(S.fluid1), cr(r.mask)}, ...
    {COL.expert, COL.ours}, struct('lineWidth',2)));
title({'5. result','green = doctor, red = ours'},'FontSize',11);
save_(f,'P02_how_it_works');
fprintf('  P02 done\n');


%% =====================================================================
%  P03  THE MONEY SLIDE -- our result against the doctor's tracing
%  =====================================================================
% SELECTION RULE, stated so it can be defended: one scan from each of six
% DIFFERENT patients, and within each patient the scan carrying the largest
% expert-marked lesion.  That is a criterion fixed in advance on the ground
% truth -- it is NOT "the six scans we happened to score best on", which would
% be cherry-picking and would not survive a question about it.
subjOrder = [8 6 5 1 10 4];
pick = zeros(numel(subjOrder),2);
for q = 1:numel(subjOrder)
    rowsQ = E(E.subject == subjOrder(q), :);
    [~,bq] = max(rowsQ.nFluid1);          % largest lesion for this patient
    pick(q,:) = [rowsQ.subject(bq), rowsQ.bscan(bq)];
end
f = figure('Color','w','Position',[10 10 1500 1300]);
tl = tiledlayout(f,3,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {'Our automatic result (red) against the eye specialist''s tracing (green)', ...
     'six different patients, largest lesion for each -- one fixed setting, no training data'});

for i = 1:size(pick,1)
    Si = loadBScan(pick(i,1), pick(i,2), 1);
    ri = segmentCME(Si.img, cfg, Si.stack);
    v  = ri.roi.mask & Si.gradedMask;
    Ei = evalSegmentation(ri.mask & v, Si.fluid1, v);
    [g0,g1,gc0,gc1] = tightCrop(ri, Si.fluid1);
    cc = @(X) X(g0:g1, gc0:gc1, :);
    nexttile;
    imshow(vizOverlay(cc(Si.img), {cc(Si.fluid1), cc(ri.mask & v)}, ...
        {COL.expert, COL.ours}, struct('fillAlpha',0.18,'lineWidth',2)));
    title(sprintf('patient %d   —   overlap with the specialist: %.0f%%', ...
        i, 100*Ei.dice), 'FontSize',11);
end
save_(f,'P03_results_gallery');
fprintf('  P03 done\n');


%% =====================================================================
%  P04  Accuracy, put in a context a non-specialist can judge
%  =====================================================================
% A bare "Dice 0.43" means nothing to this audience and sounds poor.  The
% honest and far more informative comparison is against how well the two human
% specialists agree WITH EACH OTHER on the same scans -- that is the ceiling
% any automatic method is working towards.  Both numbers are measured.
IO = load(fullfile(p.tables,'interobserver.mat'));      % dAB, from run_05

[~,ordSz] = sort(E.nFluid1,'descend');
big20     = ordSz(1:20);

% FAIR COMPARISON.  The ceiling must be measured on the SAME scans as our
% score.  Two specialists also agree with each other more easily on a big
% obvious cyst than on a faint one, so quoting our larger-cyst number against
% an all-scans ceiling would flatter us.  Both are computed per subset.
oursAll     = mean(ps.diceA,'omitnan');
ceilAll     = mean(IO.dAB,'omitnan');
oursBig     = mean(ps.diceA(big20),'omitnan');
ceilBig     = mean(IO.dAB(big20),'omitnan');
ceilingDice = ceilAll;

f = figure('Color','w','Position',[40 40 1250 660]);
M = [oursAll, ceilAll; oursBig, ceilBig] * 100;
b = bar(M, 'grouped');
b(1).FaceColor = COL.bar;
b(2).FaceColor = [0.5 0.5 0.5];
set(gca,'XTickLabel',{'all 78 scans','the larger cysts'},'FontSize',13);
ylabel('overlap with the reference tracing (%)','FontSize',12);
ylim([0 85]); grid on
legend({'our automatic method','two specialists against each other'}, ...
    'Location','northoutside','Orientation','horizontal','FontSize',12);
for kk = 1:numel(b)
    text(b(kk).XEndPoints, b(kk).YEndPoints, compose('%.0f%%',b(kk).YEndPoints), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',15,'FontWeight','bold');
end
title({'How close do we get?', ...
       'grey is the ceiling -- two trained specialists tracing the very same scans'}, ...
    'FontSize',13);
save_(f,'P04_accuracy_in_context');
fprintf('  P04 done\n');


%% =====================================================================
%  P05  The clinical output -- the part neither reference paper produces
%  =====================================================================
f = figure('Color','w','Position',[20 20 1550 680]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'From a contour to something a clinician can act on');

k = find(strcmp(string(ps.grade),'severe'), 1);
if isempty(k), k = 1; end
Sc = loadBScan(ps.subject(k), ps.bscan(k), 1);
rc = segmentCME(Sc.img, cfg, Sc.stack);
vc = rc.roi.mask & Sc.gradedMask;
d0 = max(1, round(min(rc.roi.ilm))-25);
d1 = min(size(Sc.img,1), round(max(rc.roi.obm))+25);

nexttile;
imshow(vizOverlay(Sc.img(d0:d1,:), {rc.mask(d0:d1,:) & vc(d0:d1,:)}, {COL.ours}, ...
    struct('fillAlpha',0.40,'lineWidth',2)));
txt = sprintf(['  fluid area     %.2f mm^2\n' ...
               '  cysts found    %d\n' ...
               '  mean diameter  %.0f \\mum\n' ...
               '  retinal thickness %.0f \\mum\n\n' ...
               '  SEVERITY:  %s'], ...
    rc.clinical.areaMm2, rc.clinical.cystCount, rc.clinical.meanDiamUm, ...
    rc.clinical.cstUm, upper(rc.clinical.grade));
text(12, 16, txt, 'Color','y','FontSize',13,'FontName','Consolas', ...
    'VerticalAlignment','top','BackgroundColor',[0 0 0],'Margin',6);
title('the numbers are produced automatically','FontSize',12);

nexttile;
cats = ["mild","moderate","severe"];
cnt  = arrayfun(@(g) sum(strcmp(string(ps.grade),g)), cats);
bb = bar(cnt,'FaceColor','flat');
bb.CData = [0.35 0.70 0.35; 0.95 0.75 0.20; 0.85 0.25 0.25];
set(gca,'XTickLabel',cellstr(cats),'FontSize',12);
ylabel('number of scans','FontSize',12); grid on
text(1:numel(cnt), cnt, string(cnt),'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontSize',15,'FontWeight','bold');
title('severity grade assigned across the dataset','FontSize',12);
save_(f,'P05_clinical_output');
fprintf('  P05 done\n');


%% =====================================================================
%  P06  Speed -- an easy, entirely honest win
%  =====================================================================
f = figure('Color','w','Position',[40 40 1150 600]);
tvals = [mean(ps.timeSec), 1.2, 18.4, 49.0];   % ours; and the three times the
                                               % 2021 paper quotes for itself,
                                               % DRLSE and the Snake model
b = bar(tvals,'FaceColor','flat');
b.CData = [COL.bar; 0.55 0.55 0.55; 0.55 0.55 0.55; 0.55 0.55 0.55];
set(gca,'XTickLabel',{'ours','published\newlinewave method','level set\newline(DRLSE)','active contour\newline(Snake)'}, ...
    'TickLabelInterpreter','tex','FontSize',12);
ylabel('seconds to process one scan','FontSize',12);
set(gca,'YScale','log'); grid on; ylim([0.1 100]);
text(1:4, tvals, compose('%.1f s',tvals), 'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontSize',13,'FontWeight','bold');
title({'Fast enough to run on an ordinary laptop', ...
       'times for the other three methods are as published by Liu et al. (2021)'}, ...
    'FontSize',13);
save_(f,'P06_speed');
fprintf('  P06 done\n');


%% =====================================================================
%  P07  Why four sweep directions -- the one design choice worth showing
%  =====================================================================
if exist(fullfile(p.tables,'table_02_direction_ablation.csv'),'file')
    A2 = readtable(fullfile(p.tables,'table_02_direction_ablation.csv'));
    f = figure('Color','w','Position',[40 40 1250 620]);
    yyaxis left
    b = bar(A2.Recall*100, 0.6, 'FaceColor',COL.bar); grid on
    ylabel('fluid correctly found (%)','FontSize',12); ylim([0 80]);
    text(1:height(A2), A2.Recall*100, compose('%.0f%%',A2.Recall*100), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',13,'FontWeight','bold');
    yyaxis right
    plot(A2.SecPerScan,'-o','LineWidth',2.5,'MarkerSize',9,'MarkerFaceColor','w');
    ylabel('seconds per scan','FontSize',12); ylim([0 1.2]);
    set(gca,'XTick',1:height(A2),'XTickLabel',{'4 directions','6 directions','8 directions'}, ...
        'FontSize',12);
    title({'Sweeping the operator in four directions is both the most accurate and the fastest', ...
           'this reproduces the choice made in the reference paper, on our data'}, ...
        'FontSize',12);
    save_(f,'P07_why_four_directions');
    fprintf('  P07 done\n');
end


%% =====================================================================
%  P08  A still of the live tool, so the demo has a fallback
%  =====================================================================
try
    app = CMEApp();
    app.SubjectDD.Value = 5; app.onSubjectChanged();
    app.BscanSlider.Value = 24; app.onScanChanged();
    drawnow; pause(2);
    % exportgraphics cannot capture uifigure UI components, so grab the
    % window off the screen instead.  getframe on a uifigure returns the
    % rendered window including every control.
    fr = getframe(app.Fig);
    imwrite(fr.cdata, fullfile(outDir,'P08_gui.png'));
    delete(app);
    fprintf('  P08 done\n');
catch ME
    fprintf(2,'  P08 (GUI screenshot) skipped: %s\n', ME.message);
end


%% ----------------------------------------------------------------------
writeDeckNotes(outDir, oursAll, oursBig, ceilAll, ceilBig, mean(ps.timeSec), height(E));
fprintf('\nPresentation set written to %s\n', outDir);


% =========================================================================
function [r0,r1,c0,c1] = tightCrop(res, gt)
%TIGHTCROP  Row and column range that frames the retina and the lesion.
% Rows follow the detected ILM/OBM; columns follow the annotated fluid with a
% generous margin, falling back to the whole width if there is none.
[M,N] = size(res.roi.mask);
r0 = max(1, round(min(res.roi.ilm)) - 28);
r1 = min(M, round(max(res.roi.obm)) + 28);
[~, cc] = find(gt);
if isempty(cc)
    c0 = 1; c1 = N;
else
    half = max(180, round(0.8*(max(cc)-min(cc))));
    mid  = round((min(cc)+max(cc))/2);
    c0 = max(1, mid-half); c1 = min(N, mid+half);
end
end


% =========================================================================
function annotation_arrow(cxy, label)
% Draw a short arrow pointing at a feature, with a readable label.
ax = gca;
x = cxy(1); y = cxy(2);
plot(x, y, 'o', 'MarkerSize', 26, 'LineWidth', 2.5, 'Color', [1 0.85 0.1]);
text(x+40, y-42, label, 'Color',[1 0.85 0.1], 'FontSize',13, 'FontWeight','bold', ...
    'BackgroundColor',[0 0 0], 'Margin',3);
end


function writeDeckNotes(outDir, oursAll, oursBig, ceilAll, ceilBig, tsec, nScans)
%WRITEDECKNOTES  One page of speaker notes, so the numbers on the slides and
% the words said out loud cannot drift apart.
fid = fopen(fullfile(outDir,'SPEAKER_NOTES.txt'),'w');
fprintf(fid, 'SPEAKER NOTES - CME segmentation deck\n');
fprintf(fid, '=====================================\n\n');
fprintf(fid, 'P01  What we are looking at\n');
fprintf(fid, '  An OCT scan is a cross-section through the retina, like a slice\n');
fprintf(fid, '  through a cake. In edema, fluid collects in pockets - the dark\n');
fprintf(fid, '  holes. A clinician currently outlines those by hand.\n\n');
fprintf(fid, 'P02  How it works\n');
fprintf(fid, '  Five steps, no machine learning, no training data. Clean the\n');
fprintf(fid, '  speckle, locate the retina, find the dark pockets inside it,\n');
fprintf(fid, '  keep the ones that look like fluid.\n\n');
fprintf(fid, 'P03  Results\n');
fprintf(fid, '  Green is the specialist. Red is our code, fully automatic.\n');
fprintf(fid, '  Six different patients, one setting used for all of them.\n\n');
fprintf(fid, 'P04  How good is that, really?\n');
fprintf(fid, '  Over ALL %d expert-annotated scans we overlap the\n', nScans);
fprintf(fid, '  specialist %.0f%% of the time. On the larger cysts, %.0f%%.\n', 100*oursAll, 100*oursBig);
fprintf(fid, '  KEY POINT - the ceiling is not 100%%. Two trained eye\n');
fprintf(fid, '  specialists tracing the SAME scans agree with each other\n');
fprintf(fid, '  only %.0f%% of the time overall, and %.0f%% on the larger cysts.\n', 100*ceilAll, 100*ceilBig);
fprintf(fid, '  So we reach about %.0f%% of what two humans achieve overall,\n', 100*oursAll/ceilAll);
fprintf(fid, '  and about %.0f%% of it on the clinically significant lesions.\n', 100*oursBig/ceilBig);
fprintf(fid, '  Both numbers measured on the same scans - a like-for-like\n');
fprintf(fid, '  comparison, not our best subset against their worst.\n\n');
fprintf(fid, 'P05  Clinical output\n');
fprintf(fid, '  Neither reference paper goes past the outline. We convert it to\n');
fprintf(fid, '  fluid area, cyst count, size and a severity grade.\n\n');
fprintf(fid, 'P06  Speed\n');
fprintf(fid, '  %.2f s per scan on a normal laptop, no GPU.\n\n', tsec);
fprintf(fid, 'P07  Why four directions\n');
fprintf(fid, '  The operator sweeps the image like a wave. Four sweep directions\n');
fprintf(fid, '  find the most fluid AND run fastest. Our data reproduces the\n');
fprintf(fid, '  choice the original authors made.\n\n');
fprintf(fid, 'IF ASKED - honest answers, all backed by docs/FINDINGS.md\n');
fprintf(fid, '  "Did you reproduce the paper exactly?"  We implemented it and it\n');
fprintf(fid, '     did not transfer to this dataset unchanged; our modifications\n');
fprintf(fid, '     were needed. The measurements are in the report.\n');
fprintf(fid, '  "Is it overfitted?"  No - leave-one-subject-out validation, where\n');
fprintf(fid, '     settings are tuned on nine patients and tested on the tenth,\n');
fprintf(fid, '     loses only ~0.015 overlap.\n');
fprintf(fid, '  "Any machine learning?"  None anywhere.\n');
fclose(fid);
end
