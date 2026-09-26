# CME Segmentation from OCT B-Scans — Project Report

**BME 404 Medical Imaging Sessional · Section A2, Group 03**
Monon Mohammad Sadim Sami (2118027) · Nabil Faruque Rafin (2118029) · Suraiya Sujana Khan (2118041)

*New document. Nothing existing was modified. This is the short version —
`FINDINGS.md` has the full evidence, `BUILTINS.md` the full function notes.*

---

## 1. What the project does

Cystoid macular edema (CME) is fluid trapped in pockets inside the retina, a
major treatable cause of vision loss. Clinicians outline those pockets by hand,
which is slow. Published automatic methods are mostly deep learning, needing a
GPU and large labelled training sets.

We built a **fully classical, training-free MATLAB pipeline** that segments CME
from a single OCT B-scan in about half a second on an ordinary laptop, and
converts the contour into a clinical read-out (fluid area, cyst count,
diameter, severity grade).

**Reference papers** (both post-2015):
- Liu, Lou, Chen, Cai & Wang (2021), *Applied Sciences* 11(14):6480 — the
  omnidirectional wave operator. Our base pipeline.
- Lou, Chen, Han, Liu, Wang & Cai (2020), *IEEE Access* 8:53678 — the wave
  potential-energy equation and Z/Q templates the 2021 paper builds on.

**Dataset**: Duke SD-OCT DME set (Chiu et al. 2015). 10 patients × 61 B-scans =
610 images. Only **110 were ever graded** by the two experts; of those **78
contain fluid** (our evaluation set) and **19 are fluid-free per both graders**.
The other 500 have all-NaN annotations — never examined, so nothing we output
there can be called an error.

---

## 2. The pipeline — five stages

| # | Stage | What happens |
|---|---|---|
| 1 | **Preprocess** | median filter removes speckle; CLAHE fixes uneven brightness |
| 2 | **ROI** | row-mean-gray curve finds the retina; ILM/OBM traced per A-scan; search restricted to the inner retina where CME occurs |
| 3 | **Wave operator** | fluid-mechanics analogy: gray value → fluid height, so a boundary is where "potential energy" peaks. Swept in 4 directions (0, π/2, π, 3π/2) to close a contour around an elliptical cyst |
| 4 | **Integration** | the 4 directional contours are fused; cyst bodies recovered by hysteresis thresholding; components screened by darkness, size, position and surround contrast |
| 5 | **Clinical output** | area (mm²), cyst count, equivalent diameter, central subfield thickness → mild / moderate / severe |

The wave equation implemented literally from the papers:

```
wave = φv + φg           φg = Σ I(i,j)·g(i,j) / MAX      (gravitational)
                         φv = v · vq · σ                  (kinetic)
                         σ  = exp(−(I/255δ)²),  δ = median/mean in the 3×3 Z template
boundary point iff  K = Kb/Kc > 1  AND  C = Sc/Sb > 1     (correction equation)
```

---

## 3. Code map

**`src/` — all algorithm code. The GUI and every script call these; no logic is duplicated.**

| File | Purpose |
|---|---|
| `cmeConfig.m` | every tunable parameter in one struct, each tagged [PAPER] / [DUKE] / [LIT] / [OURS] |
| `loadBScan.m` | one B-scan + both expert masks + the graded window + neighbour slices |
| `buildEvalSet.m` | indexes all 610 B-scans, caches the result |
| `preprocessOCT.m` | Stage 1 — denoise + contrast |
| `extractROI.m` | Stage 2 — row-mean-gray, ILM/OBM, inner band |
| `waveMaps.m` | the wave potential-energy and correction equations (one direction) |
| `omniWaveContours.m` | the direction-adjustment function — rotates and re-runs for 4 directions |
| `integrateContours.m` | Stage 4 — fusion and component screening |
| `niblackThreshold.m` | Niblack local threshold, written from scratch (MATLAB has none) |
| `clinicalQuantify.m` / `gradeSeverity.m` | Stage 5 — measurements and severity grade |
| `evalSegmentation.m` | Dice, precision, recall, F1, accuracy |
| `segmentCME.m` | **top-level entry point** — runs all five stages |
| `runBatch.m` / `runLOSO.m` | batch evaluation and leave-one-subject-out validation |
| `vizOverlay.m`, `figStyle.m`, `figTitle.m`, `dukePaths.m` | figure and path helpers |
| `segmentCME3D.m` | optional neighbour-slice variant (tested, not adopted — see `FINDINGS_volumetric.md`) |

**`scripts/`** — `run_00` … `run_11` produce every figure and table;
`run_ALL.m` runs the lot (~38 min); `run_P_presentation.m` builds the slide
deck; `dev_*.m` are the measurement scripts behind the findings.

**`gui/CMEApp.m`** — the demo GUI. Loads any B-scan, runs the pipeline, shows
the overlay, 16 inspectable intermediate stages, the clinical read-out and the
accuracy against both graders. Contains no algorithm code.

---

## 4. MATLAB functions — what each one actually does, and why it is used here

| Function | What it computes | Why it is the right tool here |
|---|---|---|
| `medfilt2` | replaces each pixel by the **median** of its neighbourhood | a rank filter, so it *discards* speckle outliers instead of averaging them in — and a step edge survives, which the wave operator needs |
| `imgaussfilt` | convolution with a Gaussian kernel (linear low-pass) | the baseline's denoiser; averages *across* edges, which is exactly the weakness we replaced |
| `adapthisteq` | CLAHE — histogram-equalises each tile separately, clips each bin, interpolates between tiles | handles the uneven illumination of a B-scan. Tiles must be **larger than a cyst**, so 4×4 not 8×8 |
| `imadjust` + `stretchlim` | maps the 1st–99th intensity percentiles to [0,1] with a power law | percentile limits stop a few saturated speckle pixels dictating the mapping |
| `fspecial('gaussian',[3 3])` | a 3×3 Gaussian kernel summing to 1 | **is** the weighting function g(i,j) in the φg equation |
| `imfilter(...,'corr')` | correlation (kernel **not** flipped), borders replicated | our kernels are deliberately asymmetric to reach the Q template *ahead* of Z; convolution would flip them and reach behind — silently reversing the operator |
| `diff(D,1,1)` | first differences down the columns | the difference approximation to the fluid velocity in the papers |
| `rot90` | rotates by 90° by **re-indexing** — no interpolation | implements the direction-adjustment function exactly; this is why 4 directions do not alter the image at all |
| `imrotate(...,'bilinear')` | arbitrary-angle rotation, interpolated | only for the 6/8-direction ablation — it *invents* pixel values, the paper's "image authenticity" objection made concrete |
| `graythresh` | **Otsu**: picks the threshold maximising between-class variance | fully automatic tissue/background split. Caveat: it always returns *a* threshold, even for a unimodal histogram |
| `movmean` / `movmedian` | sliding-window mean / rank filter along a vector | `movmean` smooths the row-mean curve so speckle cannot create spurious crossings; `movmedian` repairs layer traces (a mean would be dragged by the outlier) |
| `interp1(...,'extrap')` | piecewise-linear interpolation and extension | bridges A-scans where no boundary was found |
| `bwlabel` / `bwconncomp` | label each connected group of foreground pixels | lets every candidate cyst be screened on its own |
| `regionprops` | area, mean intensity, bounding box, solidity of each labelled region | supplies every screening feature; `BoundingBox(4)` is the axial extent = "thickness" |
| `imfill(...,'holes')` | floods from the border, keeps what it cannot reach | turns a closed cyst wall into a solid region |
| `imclose` / `imopen` | dilate-then-erode / erode-then-dilate | bridge gaps in contour arcs; delete 1–2 px specks |
| `imreconstruct` | repeatedly dilates a marker, clipped to a mask | implements **hysteresis thresholding** — keeps only loose regions containing a strict seed, the same idea Canny uses |
| `bwperim` + `imdilate` | boundary pixels, then thickened | draws the contour overlays legibly at slide scale |
| `activecontour` (Chan–Vese) | evolves a boundary to make inside/outside each uniform | a *region* model needing no gradient — which matters because this data's edge SNR is ≈0.6. Tested; left off by default |
| `imnlmfilt`, `imbothat`, `imhmin`, `watershed` | non-local means; black top-hat; minima suppression; watershed | all tested as improvements and **rejected** with measurements — see `FINDINGS.md` §7d |

**Written by hand because MATLAB has no equivalent:**
`niblackThreshold.m` — `T = μ + k·σ` over a local window. `imbinarize(...,'adaptive')`
is *Bradley's* method (a different formula), so calling it Niblack would be wrong.
Local statistics come from box filters via `var = E[x²] − E[x]²`.

**No machine learning anywhere.** Every threshold is fixed, derived from the
image's own statistics at run time, or measured once from the annotations and
reported as such.

---

## 5. Results

Scored on all 78 fluid-positive B-scans, inside the retina and inside the
window the graders actually annotated:

| | Dice | Precision | Recall | s/scan |
|---|---|---|---|---|
| Our pipeline vs expert A | 0.426 | 0.393 | 0.604 | 0.46 |
| Our pipeline vs expert B | 0.465 | 0.434 | 0.557 | 0.46 |
| **Leave-one-subject-out (the number to quote)** | **0.411** | 0.394 | 0.581 | 0.46 |
| *Expert A vs Expert B — the ceiling* | *0.576* | — | — | — |

**Read the last row first.** Two trained specialists tracing the same scans
agree with each other only 58% of the time (74% on the larger cysts). Our 43%
is therefore about **74% of what two humans achieve**, and **82%** on the
larger, clinically significant lesions — not 43% of a perfect score.

Leave-one-subject-out re-tunes the parameters on nine patients and scores the
tenth: the drop is only **0.015 Dice**, and nine of the ten folds chose
identical parameters, so the method is not overfitted.

**Direction ablation** reproduces the paper's own conclusion on our data: 4
directions gives the best recall (0.580 vs 0.549 for 8) at half the time, and
`rot90` is exact so it never interpolates the image.

---

## 6. Honest limitations

1. **The local edge signal is weak.** Measured over 1277 true cyst-wall
   crossings, the gray step is ~10 levels against a speckle σ of ~15 — an edge
   SNR of **0.66**. Below 1, no purely local edge operator can close a contour
   reliably. This drove the whole design.
2. **The paper's size filter excludes 89% of real cysts.** Keeping only
   components 0.5–1.5× the retinal thickness admits just 10.6% of the 443
   expert-marked Duke cysts (median cyst = 0.23×).
3. **Nine improvements were tried and measured; none beat 0.426** — including
   using neighbouring slices, which fails because our false positives are
   *systematic* (the outer nuclear layer appears in every slice, so it is
   backed by neighbours as strongly as a real cyst: 94% vs 96%, AUC 0.528).
4. **Performance scales with lesion size**: Dice 0.6–0.8 on large cyst
   clusters, near zero where the grader marked a few hundred pixels.

---

## 7. How to run

```matlab
cd 'C:\path\to\BME404-CME-Segmentation'   % wherever you unzipped it
addpath('src','gui')

cmeGUI                          % the demo GUI

S   = loadBScan(5, 24);         % or one scan from the command line
res = segmentCME(S.img, cmeConfig('modified'), S.stack);
imshow(vizOverlay(S.img, {S.fluid1, res.mask}, {[0 1 0],[1 0 0]}))
disp(res.clinical.grade)

addpath('scripts'); run_ALL     % regenerate every figure and table (~38 min)
```

**Other documents**: `FINDINGS.md` (all measurements and failed attempts),
`BUILTINS.md` (full function notes), `VIVA_NOTES.md` (likely questions),
`FINDINGS_volumetric.md` (the neighbour-slice experiment).
