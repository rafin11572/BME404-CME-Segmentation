function p = dukePaths()
%DUKEPATHS  Central place that knows where everything on disk lives.
%
% OUTPUT
%   p.root, p.data, p.results, p.figures, p.tables, p.src, p.gui
%
% Every script calls this instead of hard-coding paths, so the project can be
% moved or zipped without editing a dozen files.

here   = fileparts(mfilename('fullpath'));   % ...\<root>\src
p.root    = fileparts(here);
p.src     = here;
p.data    = fullfile(p.root, 'data');
p.results = fullfile(p.root, 'results');
p.figures = fullfile(p.results, 'figures');
p.tables  = fullfile(p.results, 'tables');
p.gui     = fullfile(p.root, 'gui');
p.docs    = fullfile(p.root, 'docs');

for f = {'results','figures','tables','docs'}
    if ~exist(p.(f{1}), 'dir'), mkdir(p.(f{1})); end
end
end
