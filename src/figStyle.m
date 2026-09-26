function figStyle(f, name, alsoFig)
%FIGSTYLE  Save a figure to results/figures with a consistent look and name.
%
%   figStyle(f, 'fig_03_roi_extraction')
%
% Every deliverable figure goes through this so the whole set shares one
% resolution, background and naming convention, and so the manifest stays in
% step with what is actually on disk.
%
% INPUTS
%   f        figure handle
%   name     base filename, no extension (convention: fig_NN_shortname)
%   alsoFig  true to also save a .fig for later editing (default false)

if nargin < 3, alsoFig = false; end
p = dukePaths();
set(f, 'Color', 'w', 'InvertHardcopy', 'off');
drawnow;   % force a full render first: exporting a tiledlayout that still has
           % pending colorbar layout updates trips a listener bug in R2021b

% 200 DPI on the generous canvas sizes used here gives images of roughly
% 2000-3000 px on the long edge, which drops into a slide at full width
% without visible softening.
try
    exportgraphics(f, fullfile(p.figures, [name '.png']), 'Resolution', 200, ...
        'BackgroundColor', 'w');
catch
    % Fallback for the same R2021b colorbar/export interaction.
    print(f, fullfile(p.figures, [name '.png']), '-dpng', '-r200');
end
if alsoFig
    savefig(f, fullfile(p.figures, [name '.fig']));
end
fprintf('  saved %s.png\n', name);
end
