# Technical findings — what we learned building this pipeline

BME 404 Medical Imaging Sessional, Section A2, Group 03
*Interpretable, Low-Resource Cystoid Macular Edema Segmentation and Severity
Grading from OCT B-Scans*

This file records the measurements that drove our design decisions. Everything
here is reproducible from the scripts in `scripts/`. It exists so that the
report and the viva answers are backed by numbers we actually measured, rather
than by assertions.

---

## 1. The dataset, as it actually is on disk

| Fact | Value |
|---|---|
| Subjects | 10 (`Subject_01.mat` … `Subject_10.mat`) |
| B-scans total | 610 (61 per subject) |
| **B-scans that were ever graded** | **110** (11 per subject) |
| — of those, with expert-A fluid | **78** ← the evaluable set |
| — of those, fluid-free per *both* graders | **19** ← the only true negatives |
| **B-scans never graded at all** | **500** (`manualFluid1` is entirely NaN) |
| B-scans with expert-B fluid | 86 |
| Image size | 496 × 768, `double`, gray 0–255 |
| Axial pixel pitch | 3.87 µm |
| Lateral pixel pitch | 10.94–11.98 µm (we use 11.46) |

Two things about the annotations that are **not** obvious and that we had to
discover by inspection:

1. **Only 110 of the 610 B-scans were ever graded.** On the other 500,
   `manualFluid1` is entirely NaN — those scans were never looked at. They are
   **ungraded, not fluid-free.** This matters twice over:
   *(a)* an earlier version of our robustness check sampled 120 "fluid-free"
   B-scans at random and reported a false-positive rate against them; almost
   all of those scans had no ground truth at all, so that number was
   meaningless and has been recomputed against the 32 genuinely graded
   fluid-free scans (19 of them, once expert B is taken into account too);
   *(b)* scrolling through the volume in the GUI looks full
   of false positives because most B-scans have nothing to compare against.
   The GUI now says explicitly which of the three cases each B-scan is in.

2. **`manualFluid1` is an instance map, not a binary mask.** Its non-zero
   values run 1…15 — the graders numbered each cyst separately. So the dataset
   silently contains an expert **cyst count** as well as an expert area, and we
   use it to validate our cyst-count output (`fig_17`). Neither reference paper
   could run that check.

3. **The graders only annotated a central window of A-scans.** The extent of
   that window is given by the columns where the layer traces are finite
   (typically columns ~117–658 of 768). We verified across **all 78**
   fluid-positive B-scans that every annotated fluid pixel lies strictly inside
   that window. Scoring outside it would charge us false positives against
   ground truth that was never drawn, so all metrics are computed inside the
   graded window only.

---

## 2. The finding that shaped the whole design: the edge SNR is below 1

The wave operator is an **edge** operator — it locates boundaries from the way
the gray value ramps along the sweep direction. So we measured how strong that
ramp actually is on this data.

Averaging over **1277 cyst-floor crossings** (every column where an expert
fluid region ends, across four subjects), the profile across a cyst wall is:

| Preprocessing | step across the wall | speckle σ inside cysts | step / σ |
|---|---|---|---|
| raw | 9.9 | 14.8 | **0.66** |
| median 5×5 | 7.4 | 10.6 | 0.70 |
| median 9×9 | 5.8 | 10.5 | 0.55 |
| Gaussian σ=3 | 9.1 | 9.9 | 0.92 |
| anisotropic diffusion | 8.3 | 11.1 | 0.75 |
| median + CLAHE | 19.2 | 21.1 | 0.91 |

*(gray levels; `scripts/dev_profile.m`)*

**The local edge signal-to-noise ratio is below 1 for every preprocessing we
tried.** No purely local edge operator can close a contour reliably under that
condition — not the wave operator, not Canny, not a level set. This is a
property of the Duke DME data, not a bug in our implementation, and it is the
single most important thing we learned.

What *does* separate cysts is **regional** intensity, not the local edge:

- cyst gray 85.5 ± 29.5 vs non-cyst retina 141.0 ± 41.2
- per-pixel separability AUC **0.84–0.93** per B-scan

So the information is there; it is just not in the local gradient. That is why
our modified pipeline recovers the cyst **body** as a low-potential-energy
basin and uses the wave operator's contours as a fence and as the retinal
boundary detector, rather than trying to fill an arc.

---

## 3. The published cyst-size filter removes almost every real cyst

Liu et al. (2021) §2.4 keeps only connected components between **0.5× and 1.5×**
the average retinal thickness. We measured all **443 expert-marked cysts** in
the evaluable set:

| Statistic | cyst thickness ÷ retinal thickness |
|---|---|
| p5 | 0.093 |
| p25 | 0.161 |
| **median** | **0.230** |
| p75 | 0.356 |
| p95 | 0.674 |
| max | 1.507 |

**Only 10.6% of real cysts fall inside the paper's window. 95.3% fall inside
ours (0.08–1.2×).** The median cyst is a quarter of the retinal thickness, not
half of it. With the paper's literal value the pipeline returns an empty mask
on essentially every scan — which is exactly what our baseline run shows.

This is reported as a sensitivity analysis (`fig_08`, `table_03`) rather than
quietly patched.

---

## 4. Where cysts actually sit (the "position feature")

### 4a. Measured against the expert LAYER traces

The Duke set also ships eight expert layer traces, so we can say exactly which
retinal layer each annotated fluid pixel falls in (345,926 pixels over the 78
evaluable B-scans):

| Layer band | share of expert fluid |
|---|---|
| ILM – NFL/GCL | 0.0% |
| NFL/GCL – IPL/INL | 2.5% |
| IPL/INL – INL/OPL | 10.4% |
| INL/OPL – OPL/ONL | 0.3% |
| **OPL/ONL – ISM/ISE (the outer nuclear layer)** | **51.9%** |
| ISM/ISE – OS/RPE | 0.1% |
| OS/RPE – BM | 0.0% |

Two consequences, both of which changed the pipeline:

- **Over half of all DME fluid sits in the ONL.** Rashno et al. (2018)
  describe discarding "the regions between middle layers" on this same
  dataset; whatever they mean by it, discarding the ONL wholesale would delete
  half the ground truth, so we do not.
- **Only 0.1% of fluid lies below the ISM/ISE boundary.** So the correct outer
  limit of the search is the top of the ellipsoid zone, not a fraction of the
  retinal thickness. Switching to that boundary raises the fraction of expert
  fluid inside the search band from **91.2% to 100.0%**, for 33% more search
  area. We estimate the ISM/ISE boundary as the RPE peak minus 11 px, having
  measured that offset to be 9–14 px (10th–90th percentile) on the 110 graded
  B-scans.

### 4b. Measured in relative depth

The paper screens by "position feature" without saying what the position rule
is. We measured it. Relative depth 0 = ILM, 1 = OBM:

| p5 | p25 | median | p75 | p95 |
|---|---|---|---|---|
| 0.15 | 0.35 | 0.45 | 0.53 | 0.65 |

False candidates, by contrast, concentrate at relative depth ≈ 0.88 — the
outer nuclear layer and the region just above the RPE, which are genuinely
hypo-reflective but are not fluid. Restricting the search to the inner retina
raised precision from 0.13 to 0.47 in one step, and is the single most
effective screening rule in the pipeline.

---

## 5. CLAHE tile size must exceed the lesion size

Our first CLAHE setting used 8×8 tiles. On a 496×768 frame that is a tile of
62×96 px — **the same size as a cyst**. Local histogram equalisation then
equalises the cyst against its own interior and flattens exactly the contrast
we need. Moving to 4×4 tiles (124×192 px) fixed it. `fig_14` shows the effect
on the cyst/retina intensity separation directly.

This is a good example of a built-in function doing precisely what it is
documented to do, and that being the wrong thing for the task.

---

## 6. Two things the papers state that do not survive contact with the data

**(a) "Accuracy" in the 2021 paper is precision.** The abstract quotes
"accuracy … 88.8%", but 0.888 appears in Table 1 on the row labelled
*Precision*. True accuracy — (TP+TN)/total — is ~93% for us simply because
most of an OCT frame is empty vitreous. We report both and label them, and we
compute TN inside the retina only so the number is not meaningless.

**(b) Dice and F1 are the same quantity.** For binary masks both equal
2TP/(2TP+FP+FN). The paper reports 0.811 and 0.813 as if they were different
measures; the gap comes from averaging per image versus pooling over pixels.
We report which convention is used instead of presenting two numbers.

---

## 7. The realistic ceiling: the two experts do not agree with each other

Neither paper reports inter-observer agreement. We computed it: expert A vs
expert B over the 78 evaluable B-scans (`fig_12`). Their mutual Dice is the
practical upper bound for any method scored against a single grader, and on
the scans where the graders disagree strongly our method is penalised no
matter what it outputs.

---

## 7b. There is no absolute intensity guard against false positives

We checked whether a fluid-free B-scan could be recognised simply by *not*
containing dark enough pixels. It cannot. Comparing 40 fluid-positive against
40 fluid-free B-scans, the dark tail of the inner-retina intensity
distribution is essentially identical:

| statistic (inner band) | with fluid | no fluid |
|---|---|---|
| 2nd percentile gray | 42.5 | 41.0 |
| 5th percentile gray | 47.0 | 47.0 |
| mean gray | 85.4 | 88.3 |

A healthy retina contains pixels just as dark as cyst fluid (the outer nuclear
layer and the shadows under vessels), so any threshold placed low enough to
find real fluid will also fire on a healthy scan. That is why our pipeline has
a non-zero false-positive floor on fluid-free B-scans (quantified in
`fig_18`), and it is a property of single-B-scan intensity analysis rather
than of our particular threshold.

We did find one shape rule that helps: cysts are rounded (median width/height
1.6) while dark-layer fragments are long strips (median 2.4, p90 5.5). Capping
the width/height ratio at 4 halves the false-positive area on fluid-free scans
and raises the "essentially clean" rate from 2% to 12% — but it costs 17% of
Dice (0.445 → 0.370), because real cysts in advanced DME are often wide and
flat too. We therefore leave `cfg.integrate.maxAspect` at `Inf` by default and
expose it as an explicit precision/recall trade-off rather than silently
optimising one metric at the expense of the headline one.

The obvious fix is volumetric context — a true cyst persists across adjacent
B-scans in the 61-scan volume, whereas a dark speckle region does not. That is
the clearest piece of future work and is out of scope for a 2-D pipeline.

---

## 7c. Of the clinical outputs, area is usable and cyst COUNT is not

Because `manualFluid1` numbers each cyst separately, we could validate the
clinical layer's two main outputs against the graders directly (`fig_17`):

- **Fluid area agrees** (Pearson r reported on the figure), with a systematic
  tendency for us to over-segment relative to expert A — consistent with our
  recall (0.58) exceeding our precision (0.41).
- **Cyst count does not agree at all** (r = −0.05). This is not a bug: a
  connected-component labeller merges anything that touches, whereas a human
  grader splits a confluent cyst cluster by eye into separate cysts. The two
  are simply not the same operation.

So the severity grade should be driven by area (and by CST), not by count. The
count is still displayed in the GUI and stored in the tables, but it is
labelled as unreliable and is not used for grading.

---

## 7d. Things we tried that did NOT work

Recorded so they are not re-tried, and because the negative results are
themselves informative about the data.

| Idea | Why it seemed right | What happened |
|---|---|---|
| **Morphological bottom-hat** (`imbothat`) as the cyst feature | finds pixels dark *relative to their own local surroundings*, which should cancel the retina's layered brightness and any lateral drift | **Worse.** Cyst-vs-retina AUC fell from 0.805 (plain intensity) to 0.75 at the best structuring-element size, and the best achievable Dice from 0.350 to 0.275. Cysts are *absolutely* dark, not merely locally dark, so removing the background level throws away real signal. |
| **Depth normalisation** (flatten on the ILM, subtract the median A-scan profile) | removes the layer structure so only the anomaly is left | **Worse**: mean AUC 0.871 → 0.778, best Dice 0.386 → 0.301. In advanced DME the cysts themselves distort the layers, so the median profile is contaminated by the very thing it is meant to model. |
| **Watershed on the wave relief** | the paper's own seawater/basin metaphor made flesh | Dice 0.20 at best, against 0.32 for a plain threshold at the time. Over-segments badly even with h-minima suppression. |
| **Directional-enclosure screening** (require a wave contour on all four sides of a candidate) | a true cyst is a closed curve, so it must be walled in every direction | Killed almost everything: the directional contours are ~1% dense, so hardly any real cyst has support on all four sides. |
| **Elongation cap** (`maxAspect`) | cysts are rounded (median w/h 1.42), dark layer fragments are long strips (median 3.79) | **Worse**: Dice 0.426 → 0.377 at w/h ≤ 4, → 0.281 at w/h ≤ 3. The feature really does separate components, but rejecting a component also throws away the real fluid merged inside it. |
| **Volumetric median** across B-scans k−1, k, k+1 | speckle is uncorrelated between B-scans, the retina is not, so this denoises with *no* in-plane blurring | **Worse**: Dice 0.426 → 0.341, even though it raises the physical edge SNR from 0.43 to 0.61 and the cyst/retina AUC from 0.876 to 0.897. The ground truth is traced on a *single slice*, so anything that mixes information across neighbours moves the mask away from that slice's tracing. A better image is not the same as better agreement with a per-slice annotation. |
| **Non-local means** denoising | averages wall-patch with wall-patch, so it should out-denoise a median filter without softening the boundary | **Worse**: Dice 0.426 → 0.259, and about 4× slower. |
| **Layer-aware search band** (ILM → ISM/ISE instead of a depth fraction) | anatomically correct, and it captures 100% of the expert fluid against 91.2% | **Worse**: Dice 0.426 → 0.258. The band is 33% larger and the extra territory is the outer nuclear layer, the worst false-positive region, so the recall gained is swamped by the precision lost. Anatomically correct ≠ empirically better. |
| **Vessel-shadow compensation** (normalise each A-scan by its RPE peak) | vessels cast vertical dark stripes through the retina that look like cysts | Essentially neutral (0.256 → 0.241 on the layer band). |
| **Marker-controlled watershed splitting** (Girish et al.'s approach on this dataset) | the measured failure is that components *merge* a cyst with adjacent dark tissue, so split rather than reject | **No gain**: 0.426 → 0.424 at best, and the darkness screen on the split pieces changed nothing at all across every setting — the pieces are all equally dark, so cyst and neighbour are not separable by intensity even after splitting. |

---

## 7e. Leave-one-subject-out: the pipeline is NOT overfitted

Every threshold here was chosen by looking at Duke B-scans, so the score on
those same scans is optimistic by construction. We therefore re-tuned the free
parameters on nine subjects and scored the tenth, rotating through all ten
(`src/runLOSO.m`, `fig_20`):

| | Dice | Precision | Recall |
|---|---|---|---|
| tuned **and** tested on all 78 (optimistic) | 0.426 | 0.393 | 0.604 |
| **leave-one-subject-out (honest)** | **0.411** | 0.394 | 0.581 |

**The optimism is only 0.015 Dice**, and nine of the ten folds independently
selected exactly the same parameters. The tuning is therefore not fitting
noise — the parameters are stable across patients, which is the thing a
leave-one-out protocol exists to test. **0.411 is the number the report should
quote.**

Two qualifications that belong with that number:

- **Which mean.** 0.411 is the mean over all 78 held-out B-scans. Weighting
  each *patient* equally instead gives 0.399 ± 0.125 (range 0.160 to 0.595
  across the ten folds). Both are defensible; quote one consistently and say
  which one it is.
- **What the 0.015 does and does not cover.** It is the optimism of the final
  parameter *selection* only. The pipeline structure, the choice of which
  parameters to expose, and the ranges swept were all decided while looking at
  the whole dataset, so the true optimism is larger than 0.015. An unbiased
  estimate would need a locked-away test set, which 10 patients cannot
  provide. Treat 0.411 as an upper bound for a new patient from the same
  scanner, not a guarantee.

**A bug worth recording.** The first version of this validation swept
`roi.ezOffsetPx`, which is only read when `roi.bandMode` is `'layer'` — and
the pipeline runs in `'relative'` mode. That dimension therefore changed
nothing: the grid looked three times larger than the search really was, and
the "chosen parameters" string quoted a meaningless tie-break value. The
headline numbers were unaffected — the effective 12-combination search was
genuine, and re-running with a corrected grid reproduced 0.411 / 0.426
exactly — but `runLOSO` now **detects any inert grid dimension automatically
and warns**, because this class of mistake is otherwise invisible.

Per-fold held-out Dice ranges from 0.160 (subject 3) to 0.595 (subject 7),
which is the same lesion-size dependence seen in `fig_11` rather than an
instability in the method.

---

## 7f. Where the gap to Liu et al.'s 0.811 actually comes from

The obvious question is why the paper reports Dice 0.811 on this dataset while
we report 0.426. We can answer it quantitatively, because the paper's
evaluation protocol is under-specified and we can simply score **our own,
completely unchanged output** under each plausible reading.

**First, their numbers are not a miscalculation.** Dice = 2PR/(P+R) holds to
three decimals for all three of their graders (e.g. expert A: P 0.898,
R 0.678 → 0.7726, reported 0.772). They computed a genuine overlap metric.

**What the paper does not state:** which of the 110 graded B-scans were
scored; over what region TP/FP/FN were counted; and whether any boundary
tolerance was allowed — which matters, because the paper is explicitly about
*contour* extraction, and contour metrics conventionally allow a tolerance.

Scoring our identical masks nine different ways (`fig_21`, `table_09`):

| protocol (algorithm unchanged throughout) | Dice | Precision | Recall |
|---|---|---|---|
| all 78, strict overlap, expert A — **what we report** | 0.426 | 0.399 | 0.604 |
| all 78, lesion-local columns | 0.589 | 0.616 | 0.604 |
| largest 20 lesions only | 0.610 | 0.693 | 0.569 |
| largest 20 + lesion-local + better grader | 0.663 | 0.785 | 0.583 |
| **+ 3 px boundary tolerance** | **0.783** | **0.888** | 0.711 |
| + 5 px boundary tolerance | 0.827 | 0.934 | 0.753 |
| *Liu et al. reported* | *0.811* | *0.888* | *0.750* |

**The precision matches to three decimals (0.888) and the Dice to within
0.03.** So our segmentation output is already capable of producing their
reported numbers — the difference is dominated by the yardstick, not by the
algorithm.

Three further observations, all verifiable from the paper itself:

1. **Their third grader is their own.** The Duke set ships two tracings; their
   Table 1 has three. The acknowledgements thank a grader at Tianjin Eye
   Hospital for re-marking the data. Their score against that commissioned
   tracing (0.880) is well above their score against the two Duke graders
   (0.772, 0.782), so the headline 0.811 average is pulled up by it. Against
   the Duke graders alone they report **0.777**.
2. **Their own size filter implies a selective method.** Keeping only
   components between 0.5× and 1.5× the retinal thickness — which we measured
   admits just 10.6% of real Duke cysts (§3) — would produce exactly their
   reported signature: high precision (0.89) with moderate recall (0.70–0.75),
   i.e. finding the big cysts confidently and missing the small ones.
3. **Their figures show large, isolated cysts**, consistent with evaluation
   concentrated on the easier cases.

**What we can and cannot conclude.** We can say the reported gap is mostly a
protocol difference, because our own masks reproduce their precision exactly
under a looser but entirely conventional protocol. We *cannot* prove which
protocol they used, and we cannot rule out that their implementation is also
genuinely more selective than ours. The honest statement for the report is:
*under a strict, whole-retina, per-image, fixed-grader protocol our pipeline
scores 0.426; the published 0.811 is not measured under that protocol, and is
not directly comparable to it.*

---

## 8. Honest summary of where we ended up

- Our modified pipeline reaches **Dice 0.426** against expert A over all 78
  evaluable B-scans and **0.411 under leave-one-subject-out**, at **~0.46 s per
  B-scan** (the paper reports 1.2 s).
- The two expert graders agree with each other at only **Dice 0.58**, so 0.411
  is about **71% of the achievable ceiling**, not 41% of a perfect score.
- We tried eight further improvements (§7d) and **none of them beat this**.
  That is worth stating plainly: the pipeline sits at a genuine local optimum
  for this class of method, and the eight measured failures say more about the
  data than about the tuning.
- Performance scales strongly with lesion size: Dice 0.6–0.8 on the larger
  cyst clusters, near 0 on scans where the expert marked only a few hundred
  pixels.
- Our faithful re-implementation of the published pipeline scores far lower,
  for the two structural reasons in §2 and §3.
- For context, published methods evaluated on this same Duke DME fluid task
  report Dice in the ~0.53 (Chiu et al. 2015, the dataset authors' own
  kernel-regression + graph method) to ~0.65–0.70 range (Rashno et al. 2018,
  neutrosophic sets + graph cuts). The 0.811 reported by Liu et al. is well
  above everything else published on this dataset, and we were not able to
  reproduce it.

We think the honest framing is the strongest one: a fully classical,
training-free, 0.4 s/scan pipeline that works well on clinically significant
(large) lesions, with the failure modes measured and shown.
