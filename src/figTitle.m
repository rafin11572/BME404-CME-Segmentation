function figTitle(f, tl, str, fs)
%FIGTITLE  Put a super-title above a tiledlayout without it colliding with the
%          first row of tile titles.
%
%   figTitle(f, tl, 'My title')
%   figTitle(f, tl, 'My title', 13)
%
% WHY THIS EXISTS
%   title(tl, ...) draws the layout title in space the layout does not
%   actually reserve once the tiles contain images (imshow makes each axes
%   fill its tile), so on an image-heavy figure the super-title lands on top
%   of the first row's titles.  Shrinking the layout explicitly and drawing
%   the title into the gap is deterministic and behaves the same on every
%   figure size.
%
% INPUTS
%   f    figure handle
%   tl   tiledlayout handle
%   str  the title, char or cellstr (a cellstr gives multiple lines)
%   fs   font size (default 14, or 13 for multi-line)

if nargin < 4 || isempty(fs)
    if iscell(str) && numel(str) > 1, fs = 13; else, fs = 14; end
end

nLines = 1;
if iscell(str), nLines = numel(str); end
band = 0.060 + 0.032*(nLines-1);        % vertical space to reserve

tl.Units = 'normalized';
tl.OuterPosition = [0 0 1 1-band];

annotation(f, 'textbox', [0 1-band 1 band], ...
    'String', str, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'EdgeColor','none', 'FontWeight','bold', 'FontSize', fs, ...
    'FitBoxToText','off', 'Interpreter','tex');
end
