%RUN_09_LOSO  Leave-one-subject-out validation.
%
%   fig_20_loso.png, table_07_loso.csv, table_08_loso_folds.csv
%
% Every threshold in this pipeline was chosen by looking at Duke B-scans, so
% quoting the score on those same B-scans measures how well we fitted the test
% set.  This script re-tunes the free parameters on nine subjects and scores
% the tenth, rotating through all ten.  The gap between the two numbers is the
% optimism, and the held-out number is the one the report should quote.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

base = cmeConfig('modified');

% The free parameters.  Deliberately few: the more knobs, the more a
% leave-one-out protocol flatters itself.
%
% Every entry here MUST be a parameter the current configuration actually
% reads.  An earlier version of this script swept roi.ezOffsetPx, which is
% only used when roi.bandMode is 'layer' -- and the pipeline runs in
% 'relative' mode, so that dimension did nothing at all: the grid looked
% three times larger than the search really was, and the reported "chosen
% parameters" quoted a meaningless tie-break value.  runLOSO now detects and
% warns about any such inert dimension automatically.
grid = struct( ...
    'name',   {'integrate.looseK', 'integrate.seedK', 'integrate.minContrast'}, ...
    'values', {{0.5, 0.7, 0.9, 1.1}, {1.0, 1.3, 1.6}, {0, 10, 20}});

R = runLOSO(base, grid, true);

fprintf('\n================ LEAVE-ONE-SUBJECT-OUT ================\n');
disp(R.summary);
fprintf('Tuned AND tested on everything (optimistic): Dice %.3f  (%s)\n', ...
    R.fittedAll.Dice, R.fittedAllParams);
fprintf('Optimism (fitted - held out) = %.3f\n', R.optimism);
fprintf('Distinct configurations actually evaluated: %d\n', R.effectiveCombos);
if ~isempty(R.inertParameters)
    fprintf(2,'INERT grid parameters (they changed nothing): %s\n', ...
        strjoin(cellstr(R.inertParameters), ', '));
end
fprintf('\nSpread across the ten folds:\n'); disp(R.foldSpread);
fprintf('CAVEAT: %s\n', R.caveat);
fprintf('\nPer-fold:\n'); disp(R.perFold);

writetable(R.perScan, fullfile(p.tables,'table_07_loso.csv'));
writetable(R.perFold, fullfile(p.tables,'table_08_loso_folds.csv'));
save(fullfile(p.tables,'loso.mat'), 'R');

% ---------------------------------------------------------------------
f = figure('Position',[20 20 1500 620]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {'Leave-one-subject-out validation', ...
     'parameters re-tuned on nine subjects, scored on the tenth'});

nexttile;
subj = R.perFold.heldOutSubject;
b = bar([R.perFold.trainDice, R.perFold.heldOutDice]); grid on
set(gca,'XTickLabel', arrayfun(@(s) sprintf('S%02d',s), subj, 'UniformOutput', false));
ylabel('Dice'); ylim([0 1]);
legend({'on the nine training subjects','on the held-out subject'}, ...
    'Location','northoutside','Orientation','horizontal');
title('per fold');

nexttile;
vals = [R.fittedAll.Dice, R.summary.Dice];
bb = bar(vals,'FaceColor','flat');
bb.CData = [0.8 0.5 0.2; 0.2 0.5 0.8];
% NOTE: a newline inside an XTickLabel string is split into two SEPARATE
% labels by MATLAB, which silently mis-pairs them with the bars.  Keep
% each label on a single line.
set(gca,'XTickLabel',{'tuned on all 78 (optimistic)','leave-one-subject-out (honest)'});
ylabel('Dice'); ylim([0 0.8]); grid on
text(1:2, vals, compose('%.3f',vals), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
title(sprintf('optimism from tuning = %.3f Dice', R.optimism));
% State the spread too: ten folds is a small sample and the mean alone is
% not an honest summary.
xlabel(sprintf('held-out folds: %.3f \\pm %.3f  (range %.3f - %.3f)', ...
    R.foldSpread.meanFoldDice, R.foldSpread.stdFoldDice, ...
    R.foldSpread.minFoldDice, R.foldSpread.maxFoldDice));

figStyle(f,'fig_20_loso');
fprintf('\nLOSO done.\n');
