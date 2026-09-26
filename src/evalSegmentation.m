function E = evalSegmentation(pred, truth, validMask)
%EVALSEGMENTATION  Pixel-wise segmentation metrics, defined as in the papers.
%
%   E = evalSegmentation(pred, truth)              score the whole frame
%   E = evalSegmentation(pred, truth, validMask)   score inside validMask only
%
% INPUTS
%   pred       MxN logical, our segmentation
%   truth      MxN logical, the expert manual fluid mask
%   validMask  MxN logical (optional), region over which to score
%
% OUTPUT (struct)
%   E.TP E.FP E.TN E.FN     raw pixel counts
%   E.dice                  2TP / (2TP + FP + FN)          Liu 2021 Eq. (11)
%   E.precision             TP / (TP + FP)                 Liu 2021 Eq. (12)
%   E.recall                TP / (TP + FN)                 Liu 2021 Eq. (13)
%   E.f1                    2PR/(P+R)                      Liu 2021 Eq. (14)
%   E.accuracy              (TP+TN)/(TP+TN+FP+FN)          textbook accuracy
%   E.jaccard               TP / (TP+FP+FN)
%
% =====================================================================
% TWO THINGS TO BE CAREFUL ABOUT WHEN COMPARING WITH THE PAPER
% =====================================================================
% 1. The 2021 paper's ABSTRACT says "accuracy ... 88.8%", but the value 0.888
%    appears in its Table 1 on the row labelled PRECISION.  Their "accuracy"
%    is therefore the positive predictive value, not (TP+TN)/total.  We
%    compute and report BOTH, and the comparison tables say which is which.
%    This matters: true accuracy on an OCT frame is ~99% simply because most
%    pixels are empty vitreous, so quoting it would flatter us enormously.
%
% 2. Dice and F1 are algebraically the SAME quantity for binary masks
%    (both equal 2TP/(2TP+FP+FN)).  The paper reports 0.811 and 0.813,
%    which differ only because one was averaged per-image and the other
%    pooled over all pixels.  We report per-image and pooled separately
%    rather than pretending they are two different measures.

pred  = logical(pred);
truth = logical(truth);
if nargin < 3 || isempty(validMask)
    validMask = true(size(pred));
end
validMask = logical(validMask);

p = pred(validMask);
t = truth(validMask);

E.TP = nnz( p &  t);
E.FP = nnz( p & ~t);
E.TN = nnz(~p & ~t);
E.FN = nnz(~p &  t);

den = @(x) max(x, eps);

E.dice      = 2*E.TP / den(2*E.TP + E.FP + E.FN);
E.precision =   E.TP / den(  E.TP + E.FP);
E.recall    =   E.TP / den(  E.TP + E.FN);
E.f1        = 2*E.precision*E.recall / den(E.precision + E.recall);
E.accuracy  = (E.TP + E.TN) / den(E.TP + E.TN + E.FP + E.FN);
E.jaccard   =   E.TP / den(  E.TP + E.FP + E.FN);

% Guard the degenerate case: if the ground truth is empty, recall/dice are
% undefined.  Report NaN rather than a misleading 0 or 1.
if (E.TP + E.FN) == 0
    E.recall = NaN; E.dice = NaN; E.f1 = NaN; E.jaccard = NaN;
end
if (E.TP + E.FP) == 0 && (E.TP + E.FN) == 0
    E.precision = NaN;
end
end
