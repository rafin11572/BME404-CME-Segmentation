# Every MATLAB built-in used in this project — what it does and why

Course requirement: *"If using built-in MATLAB functions, make sure you
understand their underlying working — you may be asked to explain them."*

This is the answer sheet. Each entry says **what the function actually
computes**, not just what it is used for. Grouped by pipeline stage.

---

## Stage 1 — denoising and contrast

### `medfilt2(I, [5 5], 'symmetric')`
Slides a 5×5 window over the image; the output pixel is the **median** of the
25 values in that window. Internally it sorts (or partially selects) the window
and takes the 13th value.

*Why here:* it is a **rank (order-statistic)** filter, not a linear one. A
speckle spike is an outlier, and the median discards outliers completely
instead of averaging them in. A step edge also survives: the median of a window
straddling an edge is still one of the two plateau values. That edge
preservation is exactly what the wave operator needs.
`'symmetric'` mirrors the image at the border so the frame edges are not
darkened by zero-padding.

### `median(stack, 3)` — the volumetric pre-filter
Per-pixel median across the three adjacent B-scans k−1, k, k+1.

*Why:* speckle is an interference effect and is **uncorrelated between
B-scans**, while the retina is highly correlated between them. So a median
along the third axis suppresses speckle with **no in-plane blurring at all** —
something no 2-D filter can offer. Measured effect on the quantity that limits
this pipeline (gray step across a cyst wall ÷ speckle σ): 0.43 → 0.61, and the
cyst-vs-retina AUC 0.876 → 0.897. The median rather than the mean, so a cyst
present in only one slice is not smeared into its neighbours.

### `imnlmfilt(I, 'DegreeOfSmoothing', h, ...)` — non-local means
Instead of averaging a pixel with its spatial neighbours, it averages it with
pixels **elsewhere** in a search window whose surrounding **patch** looks
similar. `DegreeOfSmoothing` sets how tolerant that patch comparison is;
`SearchWindowSize` and `ComparisonWindowSize` set where it looks and how big
the compared patches are.

*Why:* a cyst wall has many similar-looking wall patches along its length, so
the filter averages wall with wall and interior with interior. It therefore
suppresses speckle much harder than a median filter while leaving the boundary
step intact. It is also the most expensive step in the pipeline, which is why
the search window is kept modest.

### `imgaussfilt(I, sigma)` *(baseline mode)*
Convolves with a 2-D Gaussian kernel — a **linear** low-pass filter; every
output pixel is a weighted average of its neighbourhood with weights falling
off as exp(−r²/2σ²).

*Why:* good at zero-mean additive noise, but because it averages **across**
edges it blurs them. That weakness is exactly what our median-filter
modification targets.

### `adapthisteq(I, 'NumTiles', [4 4], 'ClipLimit', 0.01, 'Distribution', 'rayleigh')`
Contrast-Limited Adaptive Histogram Equalisation. It (1) splits the image into
`NumTiles` rectangular tiles, (2) histogram-equalises **each tile separately**
— that is the "adaptive" part, and it is what handles the uneven illumination
of an OCT B-scan — (3) **clips** each histogram bin at `ClipLimit` and
redistributes the clipped mass before equalising, and (4) **bilinearly
interpolates** between neighbouring tiles' mappings so no tile seams appear.

*Why the clip limit:* without it, plain adaptive histogram equalisation
massively amplifies noise in the near-uniform vitreous.
*Why `'rayleigh'`:* it shapes each tile's output histogram like a Rayleigh
distribution instead of a flat one. Rayleigh statistics are the textbook model
for fully-developed **speckle amplitude**, so it is the physically motivated
choice for OCT.
*Why 4×4 and not 8×8:* see `docs/FINDINGS.md` §5 — with 8×8 tiles a tile is
the same size as a cyst, and CLAHE then equalises away the very contrast we
need.

### `imadjust(I, stretchlim(I,[0.01 0.99]), [0 1], gamma)` *(baseline mode)*
`stretchlim` returns the 1st and 99th intensity percentiles; `imadjust` then
maps that window to [0,1] and applies a power law `out = in^gamma`.

*Why:* using percentiles rather than min/max means a handful of saturated
speckle pixels cannot dictate the mapping. `gamma < 1` is a concave curve that
lifts dark mid-tones (the retina) more than bright ones.

---

## Stage 2 — ROI extraction

### `movmean(x, k)` / `movmedian(x, k)`
Sliding-window mean / median along a vector — 1-D box filter and 1-D rank
filter respectively.

*Why:* `movmean` smooths the row-mean-gray curve so speckle cannot create a
dozen spurious crossings with the image-mean line near points A and B.
`movmedian` repairs the ILM/OBM traces: a retinal boundary is smooth, so a
column that jumped tens of pixels is a detection failure, and the median of its
neighbours is the right repair (a mean would be dragged by the outlier).

### `bwlabel(BW)` (1-D use)
Assigns every connected run of `true` a unique integer label.

*Why:* several short runs of the row-mean curve can sit above the image mean
(a bright artefact in the vitreous). The retina is by far the longest run, so
we label the runs and keep the longest.

### `graythresh(I)` — **Otsu's method**
Searches every possible threshold and returns the one that **maximises the
between-class variance** of the two resulting intensity groups — equivalently,
minimises the within-class variance. Fully automatic, no parameter.

*Why:* gives a tissue/background split without tuning. Caveat we state
explicitly: Otsu always returns *a* threshold, even for a unimodal histogram,
so it never tells you "there are no two classes here".

### `interp1(x, v, xq, 'linear', 'extrap')`
Piecewise-linear interpolation, with `'extrap'` extending the line beyond the
data.

*Why:* bridges columns where no boundary was detected and extends the trace to
the image edges, so every A-scan ends up with an ILM and OBM value.

---

## Stage 3 — the wave operator

### `fspecial('gaussian', [3 3], sigma)`
Builds a 3×3 Gaussian kernel normalised to sum to 1.

*Why:* this **is** the Gaussian weighting function *g(i,j)* of Eq. (6). Because
it sums to 1, `imfilter(D, g)` is a Gaussian-weighted local mean — a
noise-robust stand-in for the raw gray value in the gravitational potential
energy term.

### `imfilter(A, h, 'replicate', 'corr')`
Correlation (not convolution — the kernel is **not** flipped) of `A` with `h`,
with the border extended by replicating edge pixels.

*Why `'corr'`:* our kernels are deliberately asymmetric (e.g. a 5×3 kernel with
ones only on its last row, to reach the Q template two rows **ahead**).
Convolution would flip the kernel and reach two rows *behind* — the opposite
direction. Getting this wrong silently reverses the operator.

### `diff(D, 1, 1)`
First-order differences along dimension 1 (down the columns).

*Why:* it is the difference approximation to the first derivative of
displacement that both papers use in place of the fluid velocity.

### `medfilt2` again, inside σ
Supplies the **median** of the 3×3 Z template at every pixel, so that
δ = median/mean can be computed as a per-pixel map instead of in a loop.

### `rot90(A, k)`
Rotates by k×90° by **re-indexing** the array — no interpolation, no new pixel
values.

*Why it matters:* this is the direction-adjustment function for θ = 0, π/2, π,
3π/2. Because it is exact, the four-direction operator does not alter the image
at all — which is precisely the "authenticity of medical images" argument the
2021 paper uses to prefer 4 directions over 6 or 8.

### `imrotate(A, theta, 'bilinear', 'crop')`
Rotates by an arbitrary angle, computing each output pixel by **bilinear
interpolation** of the 4 nearest input pixels; `'crop'` keeps the output the
same size as the input.

*Why it matters:* used only for the 6- and 8-direction ablation. It invents
gray values that were never measured — the concrete form of the paper's
objection. We pad to a square first so nothing rotates out of frame.

---

## Stage 4 — contour integration and screening

### `imclose(BW, strel('disk', r))`
**Dilation followed by erosion** with the same structuring element. Bridges
gaps narrower than the element without permanently fattening the shapes.

*Why a disk:* cyst walls curve in every direction, and the disk is the only
isotropic structuring element.

### `imfill(BW, 'holes')`
Flood-fills inward from the image border and keeps whatever the flood could
**not** reach — i.e. every region fully enclosed by foreground.

*Why:* converts a closed cyst wall into a solid cyst region.

### `imreconstruct(marker, mask)` — morphological reconstruction
Repeatedly dilates `marker` and clips it to `mask` until nothing changes. The
result is exactly the union of the connected components of `mask` that contain
at least one marker pixel.

*Why:* this is how the hysteresis (double) threshold is implemented — the same
idea Canny uses for edge tracking. A strict threshold supplies seeds that are
almost certainly fluid; a loose threshold supplies the plausible extent; only
loose regions containing a seed survive. One global cut cannot do both jobs.

### `imopen(BW, strel('disk', 2))`
Erosion followed by dilation — removes features smaller than the element while
leaving larger ones roughly intact.

*Why:* deletes 1–2 px specks without shrinking real cysts.

### `bwareaopen(BW, n)`
Removes connected components smaller than `n` pixels.

### `bwlabel(BW, 8)` / `bwconncomp(BW, 8)`
Label every 8-connected group of foreground pixels with a unique integer, so
each candidate cyst can be tested on its own. `bwconncomp` returns the same
information as index lists without materialising a full label image.

*Why 8-connectivity here but 4 elsewhere:* 8-connectivity joins diagonally
touching pixels, which is right for grouping a cyst; 4-connectivity is used
when we need a one-pixel fence to be leak-proof.

### `regionprops(L, I, 'Area','MeanIntensity','BoundingBox','Solidity',...)`
Measures each labelled region. `MeanIntensity` requires the gray image to be
passed alongside the label matrix. `BoundingBox` is `[x y width height]`, so
element 4 is the **axial** extent (our "thickness") and element 3 the lateral
extent. `Solidity` is area ÷ convex-hull area.

### `bwperim(BW)` / `imdilate(BW, se)`
`bwperim` keeps only the boundary pixels of a region — what we draw for the
contour overlays. `imdilate` thickens that one-pixel outline so it stays
visible when the figure is scaled into a slide.

### `activecontour(I, mask, n, 'Chan-Vese', 'SmoothFactor', s, 'ContractionBias', b)`
Evolves the boundary of `mask` to minimise the **Chan–Vese** energy, which
rewards a partition whose inside and outside are each as uniform as possible.
`SmoothFactor` penalises contour length; `ContractionBias` > 0 shrinks the
contour, < 0 grows it.

*Why it is the right model here:* Chan–Vese is a **region** model — it needs no
image gradient at all. That matters enormously on this dataset, where we
measured the local edge SNR at about 0.6, i.e. the gradient is the one thing
the data does not have. A cyst, meanwhile, is exactly what Chan–Vese assumes: a
near-uniform dark region inside near-uniform brighter tissue.

We constrain the result with `imreconstruct` afterwards so the evolution can
adjust a boundary but cannot invent a brand-new cyst.

### `imbothat(I, SE)` — black top-hat *(tested and rejected)*
`imclose(I,SE) - I`. Closing with a structuring element larger than a cyst
fills the cyst in, so the difference is large wherever the image is dark
**relative to its own local surroundings**.

*Why we do not use it:* it made things worse — cyst-vs-retina AUC fell from
0.805 to 0.75 and the best achievable Dice from 0.350 to 0.275. Cysts are
*absolutely* dark, not merely locally dark, so removing the background level
throws away real signal. Kept as a selectable option
(`cfg.integrate.basinFeature = 'bothat'`) with the numbers recorded.

### `imerode(BW, SE)`
Erosion: keeps a pixel only if the whole structuring element fits inside the
foreground there. Used as the optional final shave of the mask, since our
recall exceeds our precision.

### `imhmin(I, h)` *(explored, not in the final pipeline)*
Suppresses every regional minimum shallower than `h` — the standard cure for
watershed over-segmentation.

---

## Stage 5 and evaluation

### `prctile(x, p)` / `median`, `std`, `skewness`
Standard order statistics. `prctile` sorts and interpolates; used for the
hysteresis thresholds under the `'percentile'` rule and for every distribution
quoted in `FINDINGS.md`.

### `regionprops(..., 'MajorAxisLength','MinorAxisLength')`
Fits the ellipse with the **same second moments** as the region — the standard
way to describe an elongated cyst.

**Anisotropy warning we handle explicitly:** Duke pixels are 3.87 µm axially
but ~11.46 µm laterally, a 3:1 ratio. An equivalent-circle diameter computed in
*pixels* would be meaningless, so we convert the area to mm² first and derive
the diameter from the **physical** area, d = 2√(A/π).

---

## Things we implemented ourselves because MATLAB has no built-in

### `niblackThreshold.m` — Niblack's local threshold
`T(x,y) = μ(x,y) + k·σ(x,y)` over a local window.
MATLAB has **no** Niblack. `imbinarize(...,'adaptive')` is **Bradley's** method
(local mean with a sensitivity offset) — a different formula — so using it and
calling it Niblack would be wrong. We compute the local statistics with box
filters via the computational formula `var = E[x²] − E[x]²`, which is O(1) per
pixel, and `max(...,0)` guards against tiny negative variances from
floating-point round-off in flat regions.

Note on the sign of `k`: Niblack's original text uses `k = −0.2` for extracting
**dark** text from a **light** page. Our feature of interest is **bright** on a
darker map, so we use `k > 0`, which pushes the threshold **above** the local
mean.

### `directionalEnclosure.m`
Scores how completely a candidate region is walled in, per operating direction,
by rotating the problem with `rot90` so all four tests share one code path.
