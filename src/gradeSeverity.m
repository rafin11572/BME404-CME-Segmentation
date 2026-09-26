function [grade, basis, detail] = gradeSeverity(C, cfg)
%GRADESEVERITY  Map quantitative CME measurements to mild / moderate / severe.
%
%   [grade, basis, detail] = gradeSeverity(C, cfg)
%
% INPUTS
%   C    the measurement struct built by clinicalQuantify
%   cfg  config struct (cfg.clinical.gradeMode, .areaCutMm2, .cstCutUm)
%
% OUTPUTS
%   grade   "none" | "mild" | "moderate" | "severe"
%   basis   which rule produced it
%   detail  struct with BOTH gradings, so the GUI and the report can show the
%           area-based and the thickness-based answer side by side
%
% =====================================================================
% WHERE THESE CUT-OFFS COME FROM -- READ THIS BEFORE QUOTING THEM
% =====================================================================
%
% (A) CST-BASED GRADE  --  literature anchored  [LIT]
%     Central subfield thickness is the quantity that commercial OCT software
%     reports and that clinical trials use.  The DRCR Retina Network defines
%     centre-involved DME on Spectralis at a CST of about 305 um (men) /
%     290 um (women), and uses ~400 um as the threshold above which
%     anti-VEGF rather than laser is indicated.
%       Reference: Wells JA et al., "Aflibercept, Bevacizumab, or Ranibizumab
%       for Diabetic Macular Edema" (DRCR.net Protocol T), NEJM 2015;
%       372:1193-1203, and the DRCR Protocol T subgroup analysis by baseline
%       visual acuity / CST.
%     We therefore use  CST < 250 um  = mild,  250-400 um = moderate,
%     CST >= 400 um = severe.  The lower bound is rounded down from the
%     Spectralis normative central subfield (~270 +- 20 um) so that a normal
%     retina cannot be graded moderate.
%
% (B) AREA-BASED GRADE  --  project-defined, data-driven  [OURS]
%     There is NO consensus clinical cut-off for intraretinal cyst AREA on a
%     single B-scan.  Grading schemes in the literature (e.g. Otani/Panozzo
%     style morphological grading) are qualitative -- they describe pattern
%     (sponge-like / cystoid / serous detachment), not an area in mm^2.
%     We therefore do NOT pretend to cite a source for an area threshold.
%     Instead the defaults in cfg.clinical.areaCutMm2 = [0.05 0.20] mm^2 are
%     the empirical TERTILES of the total fluid area measured across the
%     evaluable Duke B-scans by this pipeline (see run_05_clinical.m, which
%     recomputes and prints them).  They are a reproducible, stated convention
%     for THIS dataset, not a clinical standard, and they are declared tunable.
%
% Because (A) is anchored and (B) is not, the default gradeMode is 'area' for
% the lesion-burden story the project is about, but BOTH are always computed
% and both are surfaced in the GUI and in the results tables.

aCut = cfg.clinical.areaCutMm2;
cCut = cfg.clinical.cstCutUm;

% ---- area-based ------------------------------------------------------
if C.areaPx == 0
    gArea = "none";
elseif C.areaMm2 < aCut(1)
    gArea = "mild";
elseif C.areaMm2 < aCut(2)
    gArea = "moderate";
else
    gArea = "severe";
end

% ---- CST-based -------------------------------------------------------
if C.cstUm < cCut(1)
    gCst = "mild";
elseif C.cstUm < cCut(2)
    gCst = "moderate";
else
    gCst = "severe";
end

detail.areaGrade   = gArea;
detail.cstGrade    = gCst;
detail.areaMm2     = C.areaMm2;
detail.cstUm       = C.cstUm;
detail.areaCutMm2  = aCut;
detail.cstCutUm    = cCut;

switch lower(cfg.clinical.gradeMode)
    case 'area'
        grade = gArea;  basis = sprintf('fluid area %.3f mm^2 (cuts %.2f/%.2f mm^2)', C.areaMm2, aCut(1), aCut(2));
    case 'cst'
        grade = gCst;   basis = sprintf('CST %.0f um (cuts %.0f/%.0f um)', C.cstUm, cCut(1), cCut(2));
    case 'both'
        % Take the more severe of the two -- a conservative screening rule.
        order = ["none","mild","moderate","severe"];
        [~, ia] = ismember(gArea, order);
        [~, ic] = ismember(gCst,  order);
        grade = order(max(ia,ic));
        basis = sprintf('max(area %s, CST %s)', gArea, gCst);
    otherwise
        error('gradeSeverity:mode','Unknown gradeMode "%s"', cfg.clinical.gradeMode);
end
end
