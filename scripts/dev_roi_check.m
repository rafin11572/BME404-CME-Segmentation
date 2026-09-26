%DEV_ROI_CHECK  Validate the ROI extraction against the expert layer traces.
% DEVELOPMENT / MEASUREMENT script -- not part of the deliverable pipeline.
% It produced numbers quoted in docs/FINDINGS.md; kept so they are reproducible.
% Development script (not a deliverable figure).

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

subs = [4 5 6 7];
errI = []; errO = []; thO = []; thG = [];
for s = subs
    m = matfile(fullfile(p.data, sprintf('Subject_%02d.mat', s)));
    L = m.manualLayers1;
    idx = find(squeeze(sum(sum(isfinite(L),1),2)) > 0)';
    for k = idx(1:min(4,numel(idx)))
        S = loadBScan(s,k);
        work = preprocessOCT(S.img, cfg);
        roi  = extractROI(work, cfg);
        gI = S.layers1(1,:); gB = S.layers1(8,:);
        ok = isfinite(gI) & isfinite(gB);
        errI(end+1) = mean(abs(roi.ilm(ok)-gI(ok))); %#ok<SAGROW>
        errO(end+1) = mean(abs(roi.obm(ok)-gB(ok))); %#ok<SAGROW>
        thO(end+1)  = roi.thicknessPx;               %#ok<SAGROW>
        thG(end+1)  = mean(gB(ok)-gI(ok));           %#ok<SAGROW>
        fprintf('S%02d k%02d  ILMerr %5.1f  OBMerr %5.1f  thick ours %5.1f gt %5.1f\n', ...
            s,k,errI(end),errO(end),thO(end),thG(end));
    end
end
fprintf('\nMEAN  ILMerr %.1f px   OBMerr %.1f px   thick ours %.1f gt %.1f\n', ...
    mean(errI), mean(errO), mean(thO), mean(thG));

% picture of the worst case
[~,w] = max(errO);
disp('done');
