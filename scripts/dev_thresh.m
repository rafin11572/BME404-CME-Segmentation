%DEV_THRESH  Does a LOCAL threshold on the wave map find cyst walls where a
% DEVELOPMENT / MEASUREMENT script -- not part of the deliverable pipeline.
% It produced numbers quoted in docs/FINDINGS.md; kept so they are reproducible.
%             global one cannot?  Scored at the BOUNDARY level, which isolates
%             the operator from the downstream region-forming step.
clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

devSet = [1 31; 4 22; 5 24; 6 26; 8 24; 10 21];
S = cell(size(devSet,1),1);
for i=1:numel(S), S{i} = loadBScan(devSet(i,1), devSet(i,2)); end

base = cmeConfig('modified');
base.contrast.claheNumTiles = [4 4];

cases = {};
cases{end+1} = {'alpha',   struct('threshold','alpha')};
cases{end+1} = {'otsu',    struct('threshold','otsu')};
for w = [15 25 41 61]
    for k = [0.1 0.2 0.3 0.5]
        cases{end+1} = {sprintf('niblack w%d k%.1f',w,k), ...
            struct('threshold','niblack','niblackWindow',w,'niblackK',k)}; %#ok<SAGROW>
    end
end

fprintf('%-20s %7s %8s %8s %8s\n','case','dens%','bRecall','bPrec','bF1');
for ci = 1:numel(cases)
    c = base;
    fn = fieldnames(cases{ci}{2});
    for f = 1:numel(fn), c.wave.(fn{f}) = cases{ci}{2}.(fn{f}); end

    dens=[]; br=[]; bp=[];
    for i=1:numel(S)
        work = preprocessOCT(S{i}.img, c);
        roi  = extractROI(work, c);
        dirs = omniWaveContours(work, c, roi.mask);
        F = dirs.fused;
        dens(end+1) = 100*mean(F(roi.mask)); %#ok<SAGROW>

        gtB = bwperim(S{i}.fluid1) & roi.mask;
        if nnz(gtB)==0, continue; end
        % within 3 px tolerance, in both directions
        Fd = imdilate(F,   strel('disk',3));
        Gd = imdilate(gtB, strel('disk',3));
        br(end+1) = nnz(gtB & Fd) / nnz(gtB);          %#ok<SAGROW> boundary recall
        bp(end+1) = nnz(F & Gd)   / max(nnz(F),1);     %#ok<SAGROW> boundary precision
    end
    R = mean(br); P = mean(bp);
    fprintf('%-20s %7.2f %8.3f %8.3f %8.3f\n', cases{ci}{1}, mean(dens), R, P, 2*P*R/max(P+R,eps));
end
disp('dev_thresh done');
