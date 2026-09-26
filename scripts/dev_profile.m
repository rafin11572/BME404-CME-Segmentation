%DEV_PROFILE  Average profile of each wave quantity across a true cyst floor.
% DEVELOPMENT / MEASUREMENT script -- not part of the deliverable pipeline.
% It produced numbers quoted in docs/FINDINGS.md; kept so they are reproducible.
% Row offset 0 = the last fluid pixel; positive = deeper (into the tissue).
% Set the environment variable PREP to choose the working image.
clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

cfg = cmeConfig('modified');
cfg.contrast.claheNumTiles = [4 4];
PREP = getenv('PREP'); if isempty(PREP), PREP = 'clahe'; end
fprintf('PREP = %s\n', PREP);

off = -8:12;
q = {'wave','phi_v','phi_g','v','vq','sigma','C'};
acc = zeros(numel(q)+1, numel(off)); nAcc = 0;
noiseStd = [];

for sc = [1 31; 4 22; 5 24; 6 26]'
    S = loadBScan(sc(1), sc(2));
    switch PREP
        case 'raw',    work = double(S.img);
        case 'med5',   work = medfilt2(double(S.img),[5 5],'symmetric');
        case 'med9',   work = medfilt2(double(S.img),[9 9],'symmetric');
        case 'med11',  work = medfilt2(double(S.img),[11 11],'symmetric');
        case 'gauss3', work = imgaussfilt(double(S.img),3);
        case 'diff',   work = imdiffusefilt(double(S.img)/255,'NumberOfIterations',10)*255;
        otherwise,     work = preprocessOCT(S.img, cfg);
    end
    roi = extractROI(work, cfg);
    W   = waveMaps(work, cfg, roi.mask);
    gt  = S.fluid1 & roi.mask;

    % local noise estimate: std inside the cyst interiors (should be flat there)
    interior = imerode(gt, strel('disk',4));
    if nnz(interior) > 200
        noiseStd(end+1) = std(work(interior)); %#ok<SAGROW>
    end

    [M,N] = size(gt);
    for c = 1:N
        r = find(gt(:,c));
        if isempty(r), continue; end
        er = r([diff(r)>1; true]);          % end of each fluid run = the floor
        for e = er'
            rows = e + off;
            if any(rows<1) || any(rows>M), continue; end
            for k = 1:numel(q)
                acc(k,:) = acc(k,:) + W.(q{k})(rows, c)';
            end
            acc(end,:) = acc(end,:) + work(rows, c)';
            nAcc = nAcc + 1;
        end
    end
end
acc = acc / nAcc;
fprintf('averaged over %d cyst-floor crossings\n', nAcc);
fprintf('speckle std inside cysts = %.1f gray levels\n\n', mean(noiseStd));

fprintf('%6s', 'offset'); fprintf('%7d', off); fprintf('\n');
for k = 1:numel(q)
    fprintf('%6s', q{k}); fprintf('%7.2f', acc(k,:)); fprintf('\n');
end
fprintf('%6s', 'gray'); fprintf('%7.1f', acc(end,:)); fprintf('\n');
step = acc(end,end) - min(acc(end,:));
fprintf('\nSTEP across the wall = %.1f gray levels;  step/noise = %.2f\n', step, step/mean(noiseStd));
