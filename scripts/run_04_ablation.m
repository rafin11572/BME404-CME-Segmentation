%RUN_04_ABLATION  The two ablations the project needs:
%   (A) number of operating directions, 4 vs 6 vs 8   (mirrors Liu 2021 Table 2)
%   (B) the cyst size filter, which is the paper parameter that fails hardest
%       on the Duke data -- we measured the expert cysts to show why.
%
% Outputs: table_02_direction_ablation.csv, fig_07_direction_ablation.png,
%          fig_08_thickness_sensitivity.png, table_03_thickness_sensitivity.csv

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

T = buildEvalSet(); E = T(T.evaluable,:);
scans = [E.subject, E.bscan];

%% ======================================================================
%  (A) 4 vs 6 vs 8 directions
%  ======================================================================
% The angles are spread evenly over the full circle in each case, exactly as
% the paper describes.  Note what changes physically: at 4 directions every
% rotation is a multiple of 90 degrees, so rot90 re-indexes the pixels with no
% interpolation at all.  At 6 and 8 directions the rotation must interpolate,
% which is precisely the "authenticity of medical images" objection the paper
% raises in Sec. 3.  Our implementation makes that concrete rather than
% rhetorical -- see omniWaveContours.
dirSets = {0:90:359, 0:60:359, 0:45:359};
names   = {'4 directions','6 directions','8 directions'};

rows = {};
for d = 1:numel(dirSets)
    cfg = cmeConfig('modified');
    cfg.wave.directionsDeg = dirSets{d};
    fprintf('  %s ...\n', names{d});
    t0 = tic;
    R = runBatch(cfg, scans, struct('verbose',false));
    el = toc(t0);
    s = R.summary;
    get = @(m) s.Average(strcmp(s.Metric,m));
    rows(end+1,:) = {names{d}, numel(dirSets{d}), get('Dice'), get('Precision'), ...
                     get('Recall'), get('F1-Score'), mean(R.perScan.timeSec), el}; %#ok<SAGROW>
end
Tab2 = cell2table(rows, 'VariableNames', ...
    {'Operator','nDirections','Dice','Precision','Recall','F1Score','SecPerScan','TotalSec'});
disp(Tab2);
writetable(Tab2, fullfile(p.tables,'table_02_direction_ablation.csv'));

f = figure('Position',[20 20 1400 600]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Ablation: number of operating directions (our data, all 78 evaluable B-scans)');

nexttile;
M = [Tab2.Dice Tab2.Precision Tab2.Recall Tab2.F1Score];
b = bar(M); grid on
set(gca,'XTickLabel',Tab2.Operator);
legend({'Dice','Precision','Recall','F1'},'Location','northoutside','Orientation','horizontal');
ylabel('score'); ylim([0 1]);
for k = 1:numel(b)
    xtips = b(k).XEndPoints; ytips = b(k).YEndPoints;
    text(xtips, ytips, compose('%.2f',ytips), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',8);
end
title('segmentation quality');

nexttile;
bar(Tab2.SecPerScan,'FaceColor',[0.85 0.4 0.2]); grid on
set(gca,'XTickLabel',Tab2.Operator); ylabel('seconds per B-scan');
text(1:height(Tab2), Tab2.SecPerScan, compose('%.2f s',Tab2.SecPerScan), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
title('processing time (4 directions needs no interpolation)');
figStyle(f,'fig_07_direction_ablation');

%% ======================================================================
%  (B) the cyst size filter
%  ======================================================================
% Liu 2021 Sec. 2.4 keeps only components between 0.5x and 1.5x the average
% retinal thickness.  We measured every expert-marked cyst in the Duke set to
% see whether real cysts actually live in that window.
Sst = load(fullfile(p.tables,'cystThicknessStats.mat'));   % from the measurement run
rel = Sst.rel;

lowers = [0.02 0.05 0.08 0.15 0.25 0.35 0.5];
rows = {};
for L = lowers
    cfg = cmeConfig('modified');
    cfg.integrate.thicknessRange = [L 1.5];
    R = runBatch(cfg, scans, struct('verbose',false));
    s = R.summary; get = @(m) s.Average(strcmp(s.Metric,m));
    rows(end+1,:) = {L, get('Dice'), get('Precision'), get('Recall'), ...
                     100*mean(rel >= L & rel <= 1.5)}; %#ok<SAGROW>
    fprintf('  lower bound %.2f -> Dice %.3f (%.0f%% of expert cysts survive)\n', ...
        L, rows{end,2}, rows{end,5});
end
Tab3 = cell2table(rows,'VariableNames', ...
    {'LowerBoundXRetinalThickness','Dice','Precision','Recall','PctExpertCystsKept'});
writetable(Tab3, fullfile(p.tables,'table_03_thickness_sensitivity.csv'));

f = figure('Position',[20 20 1400 600]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'Why the published cyst-size filter fails on the Duke DME data');

nexttile;
histogram(rel, 0:0.05:1.6, 'FaceColor',[0.2 0.5 0.8]); hold on
xline(0.5,'r--','LineWidth',2); xline(1.5,'r--','LineWidth',2);
xline(0.08,'g--','LineWidth',2); xline(1.2,'g--','LineWidth',2);
grid on; xlabel('cyst thickness / mean retinal thickness'); ylabel('number of expert cysts');
legend({sprintf('%d expert cysts',numel(rel)),'paper window 0.5-1.5','','our window 0.08-1.2'}, ...
    'Location','northeast');
title(sprintf('only %.0f%% of real cysts fall in the paper''s window', ...
    100*mean(rel>=0.5 & rel<=1.5)));

nexttile;
yyaxis left
plot(Tab3.LowerBoundXRetinalThickness, Tab3.Dice, '-o','LineWidth',2); hold on
plot(Tab3.LowerBoundXRetinalThickness, Tab3.Recall, '-s','LineWidth',1.5);
ylabel('score'); ylim([0 0.8]);
yyaxis right
plot(Tab3.LowerBoundXRetinalThickness, Tab3.PctExpertCystsKept,'-^','LineWidth',1.5);
ylabel('% of expert cysts inside the window');
xline(0.5,'k--'); grid on
xlabel('lower bound of the size filter (x retinal thickness)');
legend({'Dice','Recall','% expert cysts kept','paper value 0.5'},'Location','best');
title('sensitivity of the pipeline to that one parameter');
figStyle(f,'fig_08_thickness_sensitivity');

fprintf('\nAblations done.\n');
