%RUN_10_PROTOCOL  Why do we report Dice 0.43 where Liu et al. report 0.81 on
%                 the same dataset?  How much of that gap is the ALGORITHM and
%                 how much is the EVALUATION PROTOCOL?
%
%   fig_21_protocol_sensitivity.png, table_09_protocol.csv
%
% METHOD
%   Take OUR OWN, COMPLETELY UNCHANGED segmentation output and score it under
%   every evaluation protocol the paper could plausibly have used.  Nothing
%   about the algorithm changes between rows -- only the yardstick.
%
% WHY THIS IS A FAIR QUESTION TO ASK
%   Liu et al. (2021) Sec. 2.5 defines Dice, precision and recall by the
%   standard formulas, but the paper never states:
%     * which of the 110 graded B-scans were actually scored,
%     * over what region TP/FP/TN/FN were counted (whole frame? the ROI? a
%       neighbourhood of the lesion?),
%     * whether any boundary tolerance was allowed -- which matters, because
%       the paper is explicitly about CONTOUR extraction, and contour metrics
%       conventionally allow a tolerance.
%   Their reported numbers are internally consistent (Dice = 2PR/(P+R) holds
%   to three decimals for all three of their graders), so this is not a
%   miscalculation -- it is a genuine overlap metric computed under a protocol
%   we cannot read off the paper.
%
% NOTE ON THEIR THIRD GRADER
%   The Duke dataset ships two expert tracings.  The paper's Table 1 has three
%   (A, B, C); its acknowledgements thank a grader at Tianjin Eye Hospital for
%   re-marking the data.  Their score against that third grader (Dice 0.880)
%   is markedly higher than against the two Duke graders (0.772, 0.782), so
%   their headline 0.811 average is pulled up by their own commissioned
%   tracing.  Against the Duke graders alone they report 0.777.

clear; clc; close all;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();
cfg = cmeConfig('modified');

T = buildEvalSet(); E = T(T.evaluable,:);
n = height(E);

fprintf('segmenting %d B-scans once; only the scoring changes below\n', n);
P = cell(n,1); A = cell(n,1); B = cell(n,1); V = cell(n,1); sz = zeros(n,1);
for i = 1:n
    S   = loadBScan(E.subject(i), E.bscan(i), 1);
    res = segmentCME(S.img, cfg, S.stack);
    V{i} = res.roi.mask & S.gradedMask;
    P{i} = res.mask & V{i};
    A{i} = S.fluid1;  B{i} = S.fluid2;  sz(i) = nnz(S.fluid1);
end
[~, ord] = sort(sz, 'descend');
s20 = ord(1:20);

rows = {};
rows(end+1,:) = protoRow('all 78, strict overlap, expert A', P,A,B,V, 1:n, false, 'A', 0);
rows(end+1,:) = protoRow('all 78, strict overlap, expert B', P,A,B,V, 1:n, false, 'B', 0);
rows(end+1,:) = protoRow('all 78, lesion-local columns', P,A,B,V, 1:n, true, 'A', 0);
rows(end+1,:) = protoRow('largest 20 lesions only', P,A,B,V, s20, false, 'A', 0);
rows(end+1,:) = protoRow('largest 20 + lesion-local', P,A,B,V, s20, true, 'A', 0);
rows(end+1,:) = protoRow('largest 20 + lesion-local + better grader', P,A,B,V, s20, true, 'best', 0);
rows(end+1,:) = protoRow('the above + 2 px boundary tolerance', P,A,B,V, s20, true, 'best', 2);
rows(end+1,:) = protoRow('the above + 3 px boundary tolerance', P,A,B,V, s20, true, 'best', 3);
rows(end+1,:) = protoRow('the above + 5 px boundary tolerance', P,A,B,V, s20, true, 'best', 5);

Tab = cell2table(rows, 'VariableNames', {'Protocol','Dice','Precision','Recall'});
disp(Tab);
writetable(Tab, fullfile(p.tables,'table_09_protocol.csv'));

fprintf('\nLiu et al. report Dice 0.811, Precision 0.888, Recall 0.750\n');
fprintf('Our identical output, scored with a 3 px tolerance on the larger\n');
fprintf('lesions, gives Dice %.3f, Precision %.3f, Recall %.3f\n', ...
    Tab.Dice(8), Tab.Precision(8), Tab.Recall(8));

%% ----------------------------------------------------------------------
f = figure('Position',[20 20 1600 700]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
figTitle(f, tl, {'The same segmentation output, scored nine different ways', ...
     'nothing about the algorithm changes between bars -- only the yardstick'});

nexttile;
b = barh(Tab.Dice, 'FaceColor','flat');
b.CData = repmat([0.35 0.55 0.8], height(Tab), 1);
b.CData(1,:) = [0.2 0.35 0.6];          % what we report
set(gca,'YTick',1:height(Tab),'YTickLabel',Tab.Protocol,'YDir','reverse');
xlim([0 1]); grid on; xlabel('Dice');
xline(0.811,'r--','LineWidth',2);
text(0.811, 0.35, 'Liu et al. report 0.811 ','Color','r','FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','middle','Rotation',0, ...
    'BackgroundColor','w','Margin',1);
for k = 1:height(Tab)
    text(Tab.Dice(k)+0.01, k, sprintf('%.3f',Tab.Dice(k)), 'VerticalAlignment','middle','FontSize',9);
end
title('Dice under each protocol');

nexttile;
plot(Tab.Recall, Tab.Precision, 'o-','LineWidth',1.5,'MarkerFaceColor','w'); hold on
plot(0.750, 0.888, 'rp', 'MarkerSize',20,'MarkerFaceColor','r');
text(0.755, 0.888, '  Liu et al.','Color','r','FontWeight','bold');
for k = 1:height(Tab)
    text(Tab.Recall(k)+0.008, Tab.Precision(k), sprintf('%d',k), 'FontSize',9);
end
grid on; axis square; xlim([0 1]); ylim([0 1]);
xlabel('Recall'); ylabel('Precision');
title('precision / recall of the same masks');

figStyle(f,'fig_21_protocol_sensitivity');
fprintf('\nProtocol analysis done.\n');


% =========================================================================
function row = protoRow(name, P,A,B,V, sel, lesionLocal, grader, tol)
d=[]; pr=[]; rc=[];
for j = sel(:)'
    v = V{j};
    if lesionLocal
        % Score only the A-scans that actually contain expert fluid.  This is
        % the single biggest lever: it removes almost all of the territory in
        % which false positives can occur.
        v = v & repmat(any(A{j},1), size(v,1), 1);
    end
    pmask = P{j} & v;
    [dA,pA,rA] = scoreOne(pmask, A{j}, v, tol);
    [dB,pB,rB] = scoreOne(pmask, B{j}, v, tol);
    switch grader
        case 'A',    d(end+1)=dA; pr(end+1)=pA; rc(end+1)=rA;
        case 'B',    d(end+1)=dB; pr(end+1)=pB; rc(end+1)=rB;
        case 'best'
            if isnan(dA) || (~isnan(dB) && dB > dA)
                d(end+1)=dB; pr(end+1)=pB; rc(end+1)=rB;
            else
                d(end+1)=dA; pr(end+1)=pA; rc(end+1)=rA;
            end
    end
end
row = {name, mean(d,'omitnan'), mean(pr,'omitnan'), mean(rc,'omitnan')};
end


function [d,p,r] = scoreOne(pred, gt, valid, tol)
pred = pred & valid; gt = gt & valid;
if nnz(gt) == 0, d=NaN; p=NaN; r=NaN; return; end
if tol > 0
    % Boundary-tolerant scoring: a predicted pixel counts if it lies within
    % tol pixels of the truth, and a truth pixel counts if it lies within tol
    % pixels of a prediction.  Standard for contour evaluation, and much more
    % forgiving than strict region overlap.
    se = strel('disk', tol);
    gd = imdilate(gt,   se) & valid;
    pd = imdilate(pred, se) & valid;
    tp1 = nnz(pred & gd);
    tp2 = nnz(gt   & pd);
    d = (tp1 + tp2) / max(nnz(pred) + nnz(gt), eps);
    p = tp1 / max(nnz(pred), eps);
    r = tp2 / max(nnz(gt),   eps);
else
    e = evalSegmentation(pred, gt, valid);
    d = e.dice; p = e.precision; r = e.recall;
end
end
