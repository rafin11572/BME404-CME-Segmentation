%RUN_02_BATCH  Run both pipeline modes over all 78 evaluable B-scans and save
%              the metric tables (the project's Table 1 equivalent).
clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

modes = {'baseline','modified'};
R = struct();

for m = 1:numel(modes)
    fprintf('\n=== %s mode ===\n', modes{m});
    cfg = cmeConfig(modes{m});
    t0 = tic;
    out = runBatch(cfg);
    fprintf('  %d scans in %.1f s (%.2f s/scan)\n', height(out.perScan), toc(t0), mean(out.perScan.timeSec));
    disp(out.summary);
    R.(modes{m}) = out;

    writetable(out.perScan, fullfile(p.tables, sprintf('perScan_%s.csv', modes{m})));
    writetable(out.summary, fullfile(p.tables, sprintf('summary_%s.csv', modes{m})));
end

save(fullfile(p.tables,'batchResults.mat'), 'R');

% ---- combined comparison against the paper's reported numbers -----------
paper = table({'Dice';'Precision';'Recall';'F1-Score'}, [0.811;0.888;0.750;0.813], ...
    'VariableNames', {'Metric','PaperReported'});
cmp = paper;
for m = 1:numel(modes)
    s = R.(modes{m}).summary;
    v = zeros(height(cmp),1);
    for k = 1:height(cmp)
        v(k) = s.Average(strcmp(s.Metric, cmp.Metric{k}));
    end
    cmp.(modes{m}) = v;
end
cmp.Properties.VariableNames{end-1} = 'OurBaseline';
cmp.Properties.VariableNames{end}   = 'OurModified';
disp(cmp);
writetable(cmp, fullfile(p.tables,'comparison_vs_paper.csv'));

fprintf('\nmean time/scan: baseline %.2f s, modified %.2f s (paper reports 1.2 s)\n', ...
    mean(R.baseline.perScan.timeSec), mean(R.modified.perScan.timeSec));
fprintf('Saved tables to %s\n', p.tables);
