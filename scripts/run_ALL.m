function run_ALL()
%RUN_ALL  Regenerate every result and figure in the project, in order.
%
%   >> cd('<project root>'); addpath('scripts'); run_ALL
%
% Expect roughly 20-35 minutes: the direction ablation runs the whole 78-scan
% set three times and the size-filter sweep another seven times.
%
% ORDER MATTERS
%   00  sanity check         confirms the data loads and the pairing is right
%   01  measure the dataset  cyst thickness / depth statistics used later
%   02  batch evaluation     the metric tables, both pipeline modes
%   03  paper-mimicking figures
%   04  ablations            (needs 01 and 02)
%   05  our-contribution figures (needs 02)
%   06  clinical output      (needs 02)
%   07  robustness on healthy scans
%   08  explanatory schematic
%   09  leave-one-subject-out validation (the number to quote)
%   10  protocol sensitivity: why the paper reports 0.811 and we report 0.43
%   P   the slide-deck figure set (results/figures/presentation/)
%   99  write the figure manifest
%
% WHY THIS IS A FUNCTION AND NOT A SCRIPT
%   Every step script starts with `clear`.  run() evaluates a script in its
%   CALLER's workspace, so from a driver *script* the first `clear` would wipe
%   the driver's own loop variables.  Making the driver a function, and
%   running each step inside the local helper below, gives every step its own
%   scratch workspace.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'), fullfile(root,'gui'), fullfile(root,'scripts'));
here = fullfile(root,'scripts');

steps = { ...
    'run_00_sanity_check.m'
    'run_01_measure_dataset.m'
    'run_02_batch.m'
    'run_03_figures_paper.m'
    'run_04_ablation.m'
    'run_05_figures_ours.m'
    'run_06_clinical.m'
    'run_07_robustness.m'
    'run_08_schematic.m'
    'run_09_loso.m'
    'run_10_protocol.m'
    'run_P_presentation.m'
    'run_99_manifest.m' };

t0 = tic;
failed = {};
for i = 1:numel(steps)
    f = fullfile(here, steps{i});
    fprintf('\n================================================================\n');
    fprintf(' [%d/%d] %s\n', i, numel(steps), steps{i});
    fprintf('================================================================\n');
    if ~exist(f,'file')
        warning('run_ALL:missing','skipping missing %s', steps{i}); continue
    end
    try
        runOne(f);
        close all
    catch ME
        fprintf(2, 'FAILED: %s\n  %s\n', steps{i}, ME.message);
        failed{end+1} = steps{i}; %#ok<AGROW>
    end
end

fprintf('\n================================================================\n');
fprintf('All steps finished in %.1f minutes.\n', toc(t0)/60);
if isempty(failed)
    fprintf('No failures.\n');
else
    fprintf(2,'FAILED STEPS: %s\n', strjoin(failed, ', '));
end
fprintf('Figures : %s\n', fullfile(root,'results','figures'));
fprintf('Tables  : %s\n', fullfile(root,'results','tables'));
fprintf('\nLaunch the GUI with:   cmeGUI\n');
end


% =========================================================================
function runOne(f)
%RUNONE  Execute one step script in an isolated workspace.
run(f);
end
