%RUN_99_MANIFEST  Write results/figures/MANIFEST.csv describing every figure.
%
% The manifest maps each saved filename to a one-line description, so the
% figures can be dropped straight into a report or slide deck without anyone
% having to re-open them to remember what they show.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
p = dukePaths();

D = { ...
'fig_00_sanity_check.png',           'Data sanity check: raw B-scan, working image, ROI, fused contours, final mask, overlay vs expert A.'
'fig_01_pipeline_overview.png',      'Graphical abstract. All eight pipeline stages for one representative B-scan, ending with the clinical read-out panel.'
'fig_02_wave_maps.png',              'The wave operator internals: potential energy map, the seawater/wave/coastline three-region map (mirrors Fig.2 of Lou 2020), kinetic energy, the sigma regulation factor, and the boundary after the correction equation.'
'fig_03_roi_extraction.png',         'ROI extraction (mirrors Fig.3 of Liu 2021): the row-mean-gray curve against the whole-image mean, the resulting A and B cut lines, and the refined ILM/OBM band.'
'fig_04_directional_contours.png',   'Omnidirectional operator (mirrors Fig.5 of Liu 2021): the contour found in each of the four sweep directions separately, their fusion, and the final result vs expert A.'
'fig_05_segmentation_gallery.png',   'MAIN RESULT FIGURE. Final segmentation on eight B-scans spanning the range of fluid burden, with both expert tracings and per-scan Dice/precision/recall.'
'fig_06_roi_vs_expert_layers.png',   'Validation of the ROI stage: our ILM and OBM traces against the expert layer annotations, with per-scan errors in pixels.'
'fig_07_direction_ablation.png',     'Direction-count ablation, 4 vs 6 vs 8 (mirrors Table 2 / Fig.7 of Liu 2021) computed on our data: metrics and processing time.'
'fig_08_thickness_sensitivity.png',  'Why the published cyst-size filter fails here: distribution of the 443 expert cyst thicknesses against the paper window and ours, plus the pipeline sensitivity curve.'
'fig_09_baseline_vs_modified.png',   'OUR CONTRIBUTION FIGURE. Published pipeline vs our modified pipeline on the same four B-scans, side by side with the expert tracing.'
'fig_10_metrics_comparison.png',     'Grouped bar chart: the paper reported numbers, our baseline re-implementation and our modified pipeline, across Dice/Precision/Recall/F1, plus processing speed.'
'fig_11_per_scan_dice.png',          'Per-B-scan Dice for all 78 evaluable scans (sorted), and Dice against lesion size, showing that performance scales with lesion size.'
'fig_12_interobserver.png',          'Inter-observer agreement between the two expert graders, the realistic ceiling for any automatic method, with our score marked for comparison.'
'fig_13_failure_gallery.png',        'The four worst cases, shown honestly with both expert tracings, for the algorithmic-limitations discussion.'
'fig_14_preprocessing_effect.png',   'Why CLAHE tile size matters: four tile settings with the resulting cyst/retina intensity separation, showing that an 8x8 tile is cyst-sized and flattens the contrast.'
'fig_15_severity_distribution.png',  'Automatic severity grading across the 78 scans: grade counts, where the data-driven cut-offs sit, and area-based vs CST-based grading.'
'fig_16_severity_examples.png',      'Example B-scans of each grade with the computed area/count/diameter/CST printed on the image, plus the written-out grading rule and its provenance.'
'fig_17_cyst_count_validation.png',  'Our cyst count vs the graders own cyst IDs (the Duke masks number each cyst), and fluid-area agreement.'
'fig_18_false_positive_rate.png',    'Robustness on healthy tissue: false-positive area on 120 B-scans with no expert-marked fluid, against the true-positive area distribution.'
'fig_21_protocol_sensitivity.png','WHY OUR NUMBER DIFFERS FROM THE PAPER. The same unchanged segmentation output scored under nine evaluation protocols, against the 0.811 the paper reports. Shows that the gap is dominated by the evaluation yardstick rather than by the algorithm: under a lesion-local protocol with a 3 px boundary tolerance our masks give precision 0.888, matching the paper exactly.'
'fig_20_loso.png','LEAVE-ONE-SUBJECT-OUT VALIDATION. Parameters re-tuned on nine subjects and scored on the tenth, rotating through all ten, against the optimistic number obtained by tuning and testing on everything. The gap between the two is how much the tuning flattered us.'
'fig_19_wave_operator_schematic.png','Explanatory schematic of the wave operator for a non-specialist audience: the Z and Q templates drawn on a real patch of retina, the gray profile across a true cyst floor, the energy terms along the same A-scan, and the seawater/wave/coastline analogy written out.'
};

T = cell2table(D, 'VariableNames', {'filename','description'});

% keep only what actually exists, and flag anything on disk that is undocumented
exists = arrayfun(@(i) exist(fullfile(p.figures, T.filename{i}),'file')==2, (1:height(T))');
missing = T.filename(~exists);
T = T(exists,:);

onDisk = dir(fullfile(p.figures,'fig_*.png'));
undoc = setdiff({onDisk.name}, T.filename');

writetable(T, fullfile(p.figures,'MANIFEST.csv'));

fid = fopen(fullfile(p.figures,'MANIFEST.txt'),'w');
fprintf(fid, 'FIGURE MANIFEST - BME 404 A2 Group 03\n');
fprintf(fid, 'CME segmentation and severity grading from OCT B-scans\n\n');
for i = 1:height(T)
    fprintf(fid, '%-36s %s\n\n', T.filename{i}, T.description{i});
end
fclose(fid);

fprintf('Manifest written: %d figures documented.\n', height(T));

% ---- the slide-deck set, documented separately ------------------------
% These are the simplified, layman-facing versions used in the live
% presentation.  They are deliberately a SUBSET with the technical detail
% removed, so they are listed apart from the full record above.
pres = fullfile(p.figures,'presentation');
if exist(pres,'dir')
    D2 = { ...
    'P01_what_is_cme.png',        'Opener: what an OCT B-scan is, a quiet slice beside one with obvious cystoid edema, fluid pockets marked.'
    'P02_how_it_works.png',       'The pipeline in five pictures and no equations: raw scan, denoised, retina located, dark pockets picked out, final result against the specialist.'
    'P03_results_gallery.png',    'MAIN RESULT SLIDE. Our automatic contour (red) against the specialist tracing (green) for six different patients, each one that patient''''s largest lesion, with the overlap printed as a percentage.'
    'P04_accuracy_in_context.png','Accuracy against the only meaningful yardstick: how well the two expert graders agree with EACH OTHER on the very same scans. Both bars computed on the same subsets, so the comparison is like-for-like.'
    'P05_clinical_output.png',    'The clinical layer: one annotated scan showing fluid area, cyst count, diameter and severity, plus the grade distribution across the dataset.'
    'P06_speed.png',              'Processing time per scan against the three methods timed in Liu et al. (2021).'
    'P07_why_four_directions.png','Why the operator sweeps in four directions: it finds the most fluid AND runs fastest.'
    'P08_gui.png',                'A still of the live MATLAB GUI, as a fallback if the demo cannot be run.'
    };
    T2 = cell2table(D2,'VariableNames',{'filename','description'});
    keep2 = arrayfun(@(i) exist(fullfile(pres,T2.filename{i}),'file')==2, (1:height(T2))');
    T2 = T2(keep2,:);
    writetable(T2, fullfile(pres,'MANIFEST.csv'));
    fprintf('Presentation set: %d figures documented in %s\n', height(T2), pres);
end
if ~isempty(missing)
    fprintf('Not yet generated (%d): %s\n', numel(missing), strjoin(missing', ', '));
end
if ~isempty(undoc)
    fprintf('On disk but undocumented: %s\n', strjoin(undoc, ', '));
end

% ---- also list the tables ------------------------------------------------
tl = dir(fullfile(p.tables,'*.csv'));
fprintf('\n%d result tables in %s:\n', numel(tl), p.tables);
for i = 1:numel(tl), fprintf('  %s\n', tl(i).name); end
