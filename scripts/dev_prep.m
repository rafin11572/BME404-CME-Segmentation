%DEV_PREP  DEVELOPMENT / MEASUREMENT script (not part of the deliverable
%          pipeline).  Its results are recorded in docs/FINDINGS.md.
%
%          Original question: which working image actually helps the wave
%          operator?  Measures (a) contrast-to-noise between expert fluid and
%          surrounding retina, (b) fused contour density, (c) Dice.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

devSet = [1 31; 4 22; 5 24; 6 26; 8 24; 10 21];
S = cell(size(devSet,1),1);
for i=1:numel(S), S{i} = loadBScan(devSet(i,1), devSet(i,2)); end

V = {};
V{end+1} = {'raw',                @(I,c) double(I)};
V{end+1} = {'gauss1',             @(I,c) imgaussfilt(double(I),1)};
V{end+1} = {'gauss2',             @(I,c) imgaussfilt(double(I),2)};
V{end+1} = {'med5',               @(I,c) medfilt2(double(I),[5 5],'symmetric')};
V{end+1} = {'med7+g1',            @(I,c) imgaussfilt(medfilt2(double(I),[7 7],'symmetric'),1)};
V{end+1} = {'med5+g1.5',          @(I,c) imgaussfilt(medfilt2(double(I),[5 5],'symmetric'),1.5)};
V{end+1} = {'med5+clahe8',        @(I,c) clahe(medfilt2(double(I),[5 5],'symmetric'),[8 8],0.01)};
V{end+1} = {'med5+clahe4',        @(I,c) clahe(medfilt2(double(I),[5 5],'symmetric'),[4 4],0.01)};
V{end+1} = {'med5+clahe2',        @(I,c) clahe(medfilt2(double(I),[5 5],'symmetric'),[2 2],0.01)};
V{end+1} = {'med5+clahe4+g1.5',   @(I,c) imgaussfilt(clahe(medfilt2(double(I),[5 5],'symmetric'),[4 4],0.01),1.5)};
V{end+1} = {'med5+clahe2+g1.5',   @(I,c) imgaussfilt(clahe(medfilt2(double(I),[5 5],'symmetric'),[2 2],0.005),1.5)};
V{end+1} = {'gauss1+gamma',       @(I,c) gammaStretch(imgaussfilt(double(I),1),0.8)};

base = cmeConfig('modified');
base.integrate.thicknessRange = [0.1 1.5];
base.integrate.fenceDilate    = 0;

fprintf('%-20s %6s %8s %7s %7s %7s\n','variant','CNR','fused%','Dice','Prec','Rec');
for v = 1:numel(V)
    name = V{v}{1}; fn = V{v}{2};
    cnr=[]; dens=[]; dc=[]; pr=[]; rc=[];
    for i=1:numel(S)
        work = fn(S{i}.img, base);
        gray = medfilt2(double(S{i}.img),[5 5],'symmetric');
        roi  = extractROI(work, base);
        inF  = S{i}.fluid1 & roi.mask;
        outF = ~S{i}.fluid1 & roi.mask;
        if nnz(inF)>50
            cnr(end+1) = (mean(work(outF))-mean(work(inF))) / sqrt(0.5*(var(work(outF))+var(work(inF)))); %#ok<SAGROW>
        end
        dirs = omniWaveContours(work, base, roi.mask);
        dens(end+1) = 100*mean(dirs.fused(roi.mask)); %#ok<SAGROW>
        [mask,~] = integrateContours(dirs, gray, roi, base, work);
        E = evalSegmentation(mask, S{i}.fluid1, roi.mask);
        dc(end+1)=E.dice; pr(end+1)=E.precision; rc(end+1)=E.recall; %#ok<SAGROW>
    end
    fprintf('%-20s %6.2f %8.1f %7.3f %7.3f %7.3f\n', name, mean(cnr), mean(dens), ...
        mean(dc,'omitnan'), mean(pr,'omitnan'), mean(rc,'omitnan'));
end
disp('dev_prep done');

function O = clahe(I, tiles, clip)
O = 255*adapthisteq(min(max(I/255,0),1),'NumTiles',tiles,'ClipLimit',clip,'Distribution','rayleigh');
end
function O = gammaStretch(I, g)
In = min(max(I/255,0),1);
O = 255*imadjust(In, stretchlim(In,[0.01 0.99]), [0 1], g);
end
