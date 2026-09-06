# Observed joint-read multi-k interval calibration

Frozen source d2aa4d6 replayed27 completed development conditions, comprising
three independent10-Mb genomes (seeds6301–6303), three coverage levels and three
read-error tiers. All1,485 family conditions remain in the denominator. Fixed
k15/21/27/31 and a minimum20 effective reads were specified before this run.
No truth/error labels enter inference. Point estimates agree with the previous
fixed multi-k replay within floating-point tolerance. Held-out data are unused.

| Nominal coverage | Family conditions | Interval available | Contains truth | Coverage among available |
| --- | ---: | ---: | ---: | ---: |
| 1× | 495 | 9 | 9 | 100.00% |
| 5× | 495 | 66 | 60 | 90.91% |
| 20× | 495 | 424 | 401 | 94.58% |

The1× fraction describes only9 high-support cases and cannot be called perfect
low-coverage calibration. Overall499 intervals are available;920 conditions
have insufficient effective reads and66 have no finite extrapolation. At20×,
availability is85.66%, while81.01% (401/495) of all family conditions
have an available interval containing truth. Missingness is part of the result.
Families and coverage/error conditions are dependent, so these pooled fractions
are descriptive and not independent binomial observations. Seed-specific
stratified results are retained in `summary.tsv`.

The actual child completed in287.624 s/182.109 MiB peak RSS while other jobs
ran. This is not a final resource comparison. The model conditions on supplied
founders and genome size; it excludes catalogue, background-specificity,
extrapolation-bias and biological-sampling uncertainty. Approximate nominal95%
sampling intervals are not yet a production guarantee for real plant data.

## Six-panel figure legend

`figures_v2/joint_uncertainty.pdf` and `.svg` contain six panels with the exact
plotted data in `source_data.tsv`; the SVG has101 editable text nodes and no
raster image nodes. A, interval availability, sparse read support and missing
fits across495 family conditions per coverage. B, conditional truth coverage
for each independent genome; dashed line is nominal95%. C,20× conditional
coverage by unit divergence and read-error tier, with each point one genome;
dashed line is95%. D, all finite20× point estimates against planted copies;
dashed line is identity. E, interval width relative to truth against minimum
effective reads across k. F, interval missingness by planted abundance and
coverage. These panels describe the fixed conditional development model.

Both PNG versions were visually inspected. Version2 moves panel A's legend
into reserved white space without changing any source rows; the source TSVs
are byte-identical. The first version remains archived. An initial attempt to
render version2 was rejected when the automatic review service hit a usage
limit; subsequent Git reviews succeeded, a read-only diff confirmed the
layout-only change, and a retry through the same approval route succeeded.
The refusal was not bypassed. These remain development diagnostic figures;
final publication acceptance requires the broader scientific validation.
