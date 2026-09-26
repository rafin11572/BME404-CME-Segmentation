%RUN_08_SCHEMATIC  An explanatory schematic of the wave operator itself.
%
%   fig_19_wave_operator_schematic.png
%
% Built for a non-specialist audience: it shows the Z and Q templates on a real
% patch of retina, the gray-value profile the operator sees along one A-scan,
% and how the seawater / wave / coastline picture maps onto that profile.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

S = loadBScan(5,24);
res = segmentCME(S.img, cfg, S.stack);
W = res.dirs.exemplar;

% Find a column that crosses a real cyst floor, so the profile is a genuine one
gt = S.fluid1 & res.roi.mask;
colScore = sum(gt,1);
[~, col] = max(colScore);
rows = find(gt(:,col));
cystBottom = rows(end);

r0 = max(1, cystBottom-45); r1 = min(size(S.img,1), cystBottom+45);

f = figure('Position',[20 20 1650 980]);
tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, 'How the wave operator finds a cyst wall');

% ---- (a) where we are looking ------------------------------------------
nexttile;
a0 = max(1,round(min(res.roi.ilm))-25); a1 = min(size(S.img,1),round(max(res.roi.obm))+25);
imshow(vizOverlay(S.img(a0:a1,:), {S.fluid1(a0:a1,:)}, {[0 1 0]}, struct('lineWidth',2))); hold on
xline(col,'y','LineWidth',1.5);
title(sprintf('(a) the B-scan; yellow line = A-scan %d', col));

% ---- (b) the templates on a zoomed patch --------------------------------
nexttile;
c0 = max(1,col-22); c1 = min(size(S.img,2), col+22);
patch_ = res.work(r0:r1, c0:c1);
imagesc(patch_); axis image; colormap(gca, gray); hold on
cc = col - c0 + 1; rr = cystBottom - r0 + 1;
% Z template: 3x3 centred on the pixel under test
rectangle('Position',[cc-1.5, rr-1.5, 3, 3], 'EdgeColor',[0.2 0.6 1],'LineWidth',2);
% Q template: 3 wide x 2 deep, two rows AHEAD of Z
rectangle('Position',[cc-1.5, rr+1.5, 3, 2], 'EdgeColor',[1 0.6 0.1],'LineWidth',2);
plot(cc, rr, 'r+','MarkerSize',12,'LineWidth',2);
quiver(cc+8, rr-12, 0, 18, 0, 'Color','g','LineWidth',2,'MaxHeadSize',2);
text(cc+10, rr-14, 'sweep', 'Color','g','FontWeight','bold');
set(gca,'XTick',[],'YTick',[]);
title('(b) blue = Z template (3\times3), orange = Q (3\times2) ahead of it');

% ---- (c) the gray profile along that A-scan -----------------------------
nexttile;
prof = res.work(r0:r1, col);
plot(r0:r1, prof, 'k-','LineWidth',1.6); hold on
xline(cystBottom,'r--','LineWidth',2);
area(r0:cystBottom, prof(1:(cystBottom-r0+1)), 'FaceColor',[0.3 0.5 0.9],'FaceAlpha',0.25,'EdgeColor','none');
grid on; xlabel('row (depth)'); ylabel('gray value');
legend({'gray profile','cyst floor','inside the cyst'},'Location','northwest');
title('(c) what the operator actually sees');

% ---- (d) the energy terms along the same A-scan -------------------------
nexttile;
plot(r0:r1, W.phi_g(r0:r1,col), 'LineWidth',1.6); hold on
plot(r0:r1, W.phi_v(r0:r1,col), 'LineWidth',1.6);
plot(r0:r1, W.wave(r0:r1,col),  'k','LineWidth',2);
xline(cystBottom,'r--','LineWidth',2);
grid on; xlabel('row (depth)');
legend({'\phi_g  (gravitational)','\phi_v  (kinetic)','wave = \phi_v + \phi_g'},'Location','northwest');
title('(d) the wave potential energy terms');

% ---- (e) the sea analogy ------------------------------------------------
nexttile; axis off
txt = {
 '\bfTHE ANALOGY\rm'
 ''
 'The image is treated as a fluid. Each pixel is a'
 'particle; its gray value is the fluid HEIGHT h.'
 ''
 '  \phi_g = gh   gravitational potential energy'
 '               = Gaussian-weighted gray / MAX'
 ''
 '  \phi_v = v \cdot v_q \cdot \sigma   kinetic energy'
 '     v   = normalised gray difference in Z'
 '     v_q = the same in Q, just AHEAD of Z'
 '     \sigma   = exp(-(I/255\delta)^2) damps bright areas'
 ''
 'v_q is the "wind vane" of v: a lone speckle spike'
 'gives a big v but a random v_q, so the product'
 'collapses. That is how speckle is rejected'
 'without ever blurring it away.'};
text(0.02, 0.97, txt, 'Units','normalized','VerticalAlignment','top','FontSize',10.5);

% ---- (f) seawater / wave / coastline ------------------------------------
nexttile; axis off
txt = {
 '\bfSEAWATER, WAVE, COASTLINE\rm'
 ''
 'Travelling forward along the sweep direction the'
 'operator crosses three zones:'
 ''
 '  \bfseawater\rm  far from any boundary, low energy'
 '  \bfwave\rm      just before the boundary, HIGH energy'
 '  \bfcoastline\rm inside the object, low energy again'
 ''
 'The boundary is the junction between the wave'
 'and the coastline. The correction equation then'
 'vets it: a true boundary needs'
 ''
 '     K = K_b/K_c > 1   and   C = S_c/S_b > 1'
 ''
 'i.e. the slope behind is steeper than ahead, AND'
 'the tissue ahead is brighter than behind.'
 ''
 'Sweeping DOWN finds cyst FLOORS, up finds ROOFS,'
 'and sideways finds the WALLS -- which is why four'
 'directions are needed to close the contour.'};
text(0.02, 0.97, txt, 'Units','normalized','VerticalAlignment','top','FontSize',10.5);

figStyle(f,'fig_19_wave_operator_schematic');
fprintf('Schematic done.\n');
