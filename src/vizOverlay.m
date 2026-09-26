function RGB = vizOverlay(img, masks, colors, opts)
%VIZOVERLAY  Draw one or more mask outlines (and optional fills) on a B-scan.
%
%   RGB = vizOverlay(img, {m1,m2}, {[0 1 0],[1 0 0]})
%   RGB = vizOverlay(img, masks, colors, struct('fillAlpha',0.25,'lineWidth',1))
%
% INPUTS
%   img     MxN gray image on 0..255
%   masks   cell array of MxN logical masks
%   colors  cell array of 1x3 RGB triplets in [0,1]
%   opts    .fillAlpha  translucent fill strength (0 = outline only, default 0)
%           .lineWidth  outline thickness in pixels (default 1)
%           .clim       display window, default [0 255]
%
% OUTPUT
%   RGB     MxNx3 double in [0,1], ready for imshow or imwrite
%
% Written as an image-returning function rather than a plotting routine so the
% same code serves the figure scripts, the GUI and any saved PNG, and so the
% overlays look identical everywhere.

if nargin < 4, opts = struct(); end
if ~isfield(opts,'fillAlpha'), opts.fillAlpha = 0;     end
if ~isfield(opts,'lineWidth'), opts.lineWidth = 1;     end
if ~isfield(opts,'clim'),      opts.clim = [0 255];    end

I = double(img);
I = (I - opts.clim(1)) / max(opts.clim(2)-opts.clim(1), eps);
I = min(max(I,0),1);
RGB = repmat(I, 1, 1, 3);

for k = 1:numel(masks)
    m = logical(masks{k});
    if ~any(m(:)), continue; end
    c = colors{k};

    if opts.fillAlpha > 0
        for ch = 1:3
            layer = RGB(:,:,ch);
            layer(m) = (1-opts.fillAlpha)*layer(m) + opts.fillAlpha*c(ch);
            RGB(:,:,ch) = layer;
        end
    end

    % bwperim keeps only the pixels on a region's boundary, which is what we
    % want to draw; imdilate then thickens that one-pixel outline so it stays
    % visible when the figure is scaled down into a slide.
    edge = bwperim(m);
    if opts.lineWidth > 1
        edge = imdilate(edge, strel('disk', opts.lineWidth-1));
    end
    for ch = 1:3
        layer = RGB(:,:,ch);
        layer(edge) = c(ch);
        RGB(:,:,ch) = layer;
    end
end
end
