# Viva prep — likely questions and the honest answers

## "Explain the wave operator in one minute."

The image is treated as an ideal fluid. Each pixel is a particle and its gray
value is the fluid height, so the fluid potential-energy equation
phi = gh + p/rho + v^2/2 becomes an image quantity. Pressure is zero at the
free surface, so that term drops. What is left is

  wave = phi_g + phi_v

phi_g is the Gaussian-weighted gray value in a 3x3 template divided by the
image maximum: the "gravitational" term, high in bright tissue.
phi_v = v * vq * sigma is the "kinetic" term: v is the normalised gray
difference inside that 3x3 template, vq is the same quantity measured in a 3x2
template placed just AHEAD of it along the sweep, and sigma damps the whole
thing in bright regions so hyper-reflective layers are not over-segmented.

vq is the clever part. A lone speckle spike produces a large v but a random vq,
so the product collapses. Speckle is rejected without blurring it away.

The operator then sweeps in four directions (0, pi/2, pi, 3pi/2). It only fires
where the image goes dark-to-bright ALONG the sweep, so sweeping down finds
cyst floors, up finds roofs, sideways finds the walls. That is why four
directions are needed to close a contour around an elliptical cyst, and it is
the entire contribution of the 2021 paper over the 2020 one.

## "Why is your Dice 0.45 when the paper reports 0.81?"

Two measured reasons, both in docs/FINDINGS.md.

1. The local edge SNR on this dataset is below 1. We averaged 1277 crossings of
   real cyst floors: the gray step across the wall is about 10 levels against a
   speckle standard deviation of about 15. No purely local edge operator can
   close a contour under that condition.

2. The paper's cyst-size filter keeps components between 0.5x and 1.5x the
   retinal thickness. We measured all 443 expert-marked cysts: only 10.6% fall
   in that window. The median cyst is 0.23x. With the paper's literal value the
   pipeline returns an empty mask, which is exactly what our baseline shows.

For context, the dataset's own authors (Chiu et al. 2015) report about 0.53 on
this fluid task and Rashno et al. 2018 about 0.65-0.70. The 0.811 figure sits
well above everything else published on this data and we could not reproduce it.

## "What exactly did you change, and why each one?"

| Change | Reason |
|---|---|
| median filter instead of Gaussian blur | rank filter: removes speckle outliers without averaging across edges |
| CLAHE instead of gamma stretch | adapts to the uneven illumination of a B-scan; 4x4 tiles, because an 8x8 tile is cyst-sized and equalises the cyst away |
| Otsu / Niblack instead of a fixed alpha | a fixed global ratio cannot track a scan that is unusually bright or dark |
| basin + hysteresis instead of fill-the-arc | forced by finding (1) above; the contours still supply the ROI, the fence and the enclosure evidence |
| explicit inner-retina depth band | this is the paper's unstated "position feature"; we measured where cysts actually sit (p5 0.15, p95 0.65) |
| size window 0.08-1.2x instead of 0.5-1.5x | measured: keeps 95.3% of real cysts instead of 10.6% |
| clinical output layer | area, cyst count, diameter, CST, severity grade -- neither paper produces any of this |

## "Which numbers did you invent?"

None silently. Every parameter in cmeConfig carries a tag:
[PAPER] from the reference papers, [DUKE] from the dataset documentation,
[LIT] from clinical literature (cited inline), [OURS] our own tunable default.
The severity cut-offs are the clearest case: the CST bands are anchored to
DRCR.net Protocol T, and the area bands are explicitly labelled as tertiles of
our own measurements because no consensus cyst-area cut-off exists.

## "Where does it fail?"

- Small lesions. Dice is 0.6-0.8 on large cyst clusters and near zero when the
  grader marked only a few hundred pixels. fig_11 shows this directly.
- Healthy scans. There is a false-positive floor: a normal retina contains
  pixels just as dark as fluid (5th-percentile gray is 47 in both). fig_18
  quantifies it. The fix is volumetric context, which a 2-D pipeline cannot use.
- Where the two graders disagree strongly we are penalised whatever we output.
  fig_12 shows their mutual Dice, which is the realistic ceiling.

## "Is any of this trained?"

No. There is no learned or trained component anywhere. Every threshold is
either fixed, derived from the image's own statistics at run time, or measured
once from the expert annotations and reported as such.

## "Why is your GUI not an .mlapp?"

It is a programmatic uifigure class, so the entire interface is readable source
code that can be diffed and explained line by line, rather than a binary ZIP.
It contains no algorithm logic: every callback calls the same functions in src/
that the batch scripts use, so the GUI can never drift out of step with the
numbers in the report.

## Fast facts

- 610 B-scans, 78 with expert-A fluid, 86 with expert-B, 110 with layer traces
- 496 x 768, 3.87 um axial x 11.46 um lateral (3:1 anisotropic -- we convert
  area to mm^2 before deriving any diameter)
- about 0.4 s per B-scan; the paper reports 1.2 s
- 4 directions beat 6 and match 8 on Dice, with the best recall and the
  shortest time -- and rot90 is exact, so 4 directions do not interpolate the
  image at all
