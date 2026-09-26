# Using the neighbouring B-scans — what we tested and what we found

**New document. Nothing in `docs/FINDINGS.md` or any existing source file was
changed.** Everything described here lives in two new files:

| new file | what it is |
|---|---|
| `src/segmentCME3D.m` | segments k−1, k, k+1 independently, then combines the *results* |
| `scripts/run_11_volumetric.m` | the experiment, its diagnostic, `fig_22`, `table_10` |

---

## The question

An eye specialist does not grade one slice in isolation. A cyst is a
three-dimensional pocket of fluid, so it appears in the slices either side of
the one being traced. Could our pipeline use the same evidence?

## This is NOT the 3-D denoising we already rejected

`docs/FINDINGS.md` §7d records that taking a per-pixel **median across the
images** of k−1, k, k+1 before segmenting made things worse (Dice 0.426 →
0.341). That mixes the pictures, so a cyst that shifts slightly between slices
is smeared and the mask drifts away from the tracing drawn on slice k.

What is tested here is the opposite arrangement: **segment every slice
independently, never mixing the images, and use the neighbouring results as
supporting evidence.** Five combination rules were tried.

Because consecutive B-scans are ~118 µm apart azimuthally, the same cyst is
not at pixel-identical coordinates in the next slice. Neighbour masks are
therefore dilated by 3 px before comparison; without that tolerance the test
is far too strict and deletes genuine fluid.

## Result: it does not help

All 78 evaluable B-scans, scored against expert A, with a paired *t*-test on
the per-scan Dice (78 paired measurements):

| rule | Dice | Precision | Recall | vs 2-D |
|---|---|---|---|---|
| **none — the current 2-D pipeline** | **0.426** | 0.393 | 0.604 | reference |
| union (add what both neighbours agree on) | 0.421 | 0.363 | **0.685** | −0.005 (p = 0.51) |
| support ≥ 0.30 (component-level) | 0.421 | 0.396 | 0.576 | −0.005 (p = 0.15) |
| support ≥ 0.15 | 0.420 | 0.395 | 0.576 | −0.006 (p = 0.09) |
| intersect (pixel must appear in a neighbour) | 0.417 | **0.407** | 0.548 | −0.010 (p = 0.08) |
| support ≥ 0.50 | 0.417 | 0.400 | 0.564 | −0.009 (p = 0.03) |
| majority (≥ 2 of 3 slices) | 0.410 | 0.362 | 0.646 | −0.016 (p = 0.07) |
| support ≥ 0.70 | 0.410 | 0.392 | 0.554 | −0.016 (p = 0.08) |

**No rule beats the 2-D pipeline.** Every one is slightly worse, and none of
the differences is a meaningful improvement.

The rules do behave exactly as designed, which is worth noting: `union` trades
precision for recall (0.604 → 0.685) and `intersect` trades recall for
precision (0.393 → 0.407). The machinery works. It simply does not buy
accuracy.

> A caution recorded deliberately: on the single B-scan we first tried
> (subject 5, scan 24) `union` scored 0.620 against 0.553, and it looked like a
> clear win. It was not — one scan is not evidence. The full set reversed it.

## Why it fails — the decisive measurement

Neighbour agreement can only remove errors that are **random** from slice to
slice. So we measured, for every component our pipeline produced, how much
neighbour support it had, split by whether it was actually fluid:

| component type | median neighbour support |
|---|---|
| true (overlaps expert fluid) | **0.958** |
| false | **0.938** |

Separability **AUC = 0.528**, where 0.5 means no discrimination whatsoever
(112 true and 131 false components).

**Our false positives are systematic, not random.** The outer nuclear layer,
the dark bands between retinal layers and the vessel shadows all reappear in
every slice, so they are supported by the neighbours just as strongly as a
real cyst is — 94% against 96%. Volumetric consistency cannot separate the two
because the thing it tests for is true of both.

This also explains, retrospectively, why the earlier 3-D denoising failed for
a *different* reason: that attempt degraded the image, this one has a perfectly
good signal but the signal carries no information about the question.

## Recommendation

**Do not adopt it.** It costs roughly 3× the runtime (2.5 s vs 0.46 s per
B-scan, since three slices must be segmented) and it requires the whole volume
rather than a single B-scan, in exchange for no accuracy gain.

Keep it as a documented negative result. It is a good one to have: the obvious
question "did you try using the neighbouring slices, like a human does?" now
has a measured answer with a mechanism behind it, rather than a guess.

## What would actually need volumetric context

The finding above says the limitation is not slice-to-slice noise, it is that
the pipeline's errors are driven by structures that genuinely look like fluid
in every slice. Fixing that needs a feature that distinguishes a fluid pocket
from a dark layer *within* a slice — which is the same wall we hit in
`docs/FINDINGS.md` §2 (local edge SNR ≈ 0.66) and §7d (eight failed attempts).
Neighbouring slices do not supply it.
