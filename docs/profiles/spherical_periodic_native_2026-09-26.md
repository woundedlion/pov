# Periodic surface native quality experiment (2026-09-26)

Neither periodic surface qualifies for shipping at the tested 8–24 field-query
budgets. At 24 queries, 40.28% of cosine rays and 82.42% of gyroid rays remain
unresolved. These are diagnostic failures; they do not establish empty space.

The experiment uses the conservative global Lipschitz bounds and shared surface
search in [the benchmark](../../tests/periodic_surface_benchmark.cpp).
[Raw results](evidence/spherical_periodic_native_2026-09-26/results.csv) include
query counts, unresolved fractions, image differences, and native elapsed time.

## Method

Each case traces the 10,658 directions in the margin-expanded segment-0 clip of
the calibrated 288×144 display: 146 columns and 73 northern rows. The driver
chooses either half of the longitude range; this experiment measures both. The
left clip includes wrapped columns −1 through 144; the right includes columns
143 through 288. Rows are 0 through 72. These are the driver's nominal 10,368
samples plus its one-pixel scan margin. The 36 scene poses combine periods 0.7,
1.5, and 3.0; isovalues −0.5,
0, and 0.5; and four camera centers in period units: (0,0,0), (0.23,0.41,−0.17),
(0.5,0.5,0.5), and (−0.37,0.19,0.73). Each pose is measured in both clips, giving
72 windows. This includes starts in both regions and
near-tangent directions across the display. Rays start at the camera center and
end four periods away. Position tolerance is one ten-thousandth of a period.

Every field evaluation, including verification probes, counts against the total
query budget. Step and refinement caps are 1,024, so the query cap controls the
bounded experiments. The reference uses the same conservative search with 1,024
queries, rather than claiming an exact closed-form intersection oracle. Its
remaining unresolved fraction is reported explicitly.

Depth-image intensity is `1 - t / far` for verified hits and zero otherwise.
The image difference is mean absolute error in this linear unit-range intensity;
an unresolved ray rendered black remains a diagnostic failure. There is one ray
per pixel and no gradient lighting or supersampling in this native timing.

The captured executable was compiled with Clang 23.0.0git, Windows x86-64,
`-O2`, C++20, and without fast-math. Timings include trace diagnostics and depth
comparison. They are host timings, not Teensy render timings or admission data.
Build the `sdf_periodic_benchmark` native target and run its executable to repeat
the sweep.

## Measurements

Each row covers 767,376 rays. Peak evaluations equal the configured budget.

| Surface | Budget | Mean queries | Unresolved | Hit disagreement | Image MAE | Mean / peak native ms per quadrant |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cosine | 8 | 7.95 | 91.05% | 89.28% | 0.7040 | 2.346 / 5.507 |
| Cosine | 12 | 11.30 | 73.80% | 72.02% | 0.5480 | 3.331 / 6.000 |
| Cosine | 16 | 14.07 | 61.58% | 59.81% | 0.4398 | 4.105 / 5.985 |
| Cosine | 24 | 18.23 | 40.28% | 38.56% | 0.2663 | 5.266 / 7.704 |
| Cosine reference | 1,024 | 28.56 | 0.0188% | — | — | 7.814 / 11.943 |
| Gyroid | 8 | 7.30 | 89.76% | 89.03% | 0.7603 | 2.713 / 4.990 |
| Gyroid | 12 | 10.89 | 89.76% | 89.03% | 0.7603 | 3.984 / 6.200 |
| Gyroid | 16 | 14.48 | 88.61% | 87.88% | 0.7489 | 5.277 / 7.851 |
| Gyroid | 24 | 21.36 | 82.42% | 81.69% | 0.6884 | 7.624 / 10.798 |
| Gyroid reference | 1,024 | 63.22 | 0.2615% | — | — | 20.562 / 54.913 |

Cosine starts inside on 58.33% of rays; gyroid starts inside on 50.00%. An exact
zero-set start is emitted at zero depth. The reference is also bounded: grazing
and floating-point progress failures remain unresolved.

Device frame spills, firmware memory, gradient-lighting cost, and verified
supersampling cost require the separate segmented-device experiment. The native
quality result already prevents admitting these low-budget configurations.
