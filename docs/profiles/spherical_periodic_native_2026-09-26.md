# Periodic surface native quality experiment (2026-09-26)

Neither periodic surface qualifies for shipping at the tested 8–24 field-query
budgets. At 24 queries, 39.80% of cosine rays and 72.12% of gyroid rays remain
unresolved. These are diagnostic failures; they do not establish empty space.

The experiment uses the conservative global Lipschitz bounds and shared surface
search in [the benchmark](../../tests/periodic_surface_benchmark.cpp).
[Raw results](evidence/spherical_periodic_native_2026-09-26/results.csv) include
query counts, unresolved fractions, image differences, and native elapsed time.

## Method

Each case traces the 10,368 directions in one 72-column quadrant of the calibrated
288×144 display. The 36 cases combine periods 0.7, 1.5, and 3.0; isovalues −0.5,
0, and 0.5; and four camera centers in period units: (0,0,0), (0.23,0.41,−0.17),
(0.5,0.5,0.5), and (−0.37,0.19,0.73). This includes starts in both regions and
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

Each row covers 373,248 rays. Peak evaluations equal the configured budget.

| Surface | Budget | Mean queries | Unresolved | Hit disagreement | Image MAE | Mean / peak native ms per quadrant |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cosine | 8 | 7.96 | 90.26% | 89.02% | 0.7165 | 2.394 / 3.689 |
| Cosine | 12 | 11.26 | 72.22% | 70.98% | 0.5525 | 3.340 / 3.952 |
| Cosine | 16 | 13.96 | 59.90% | 58.66% | 0.4412 | 4.249 / 6.066 |
| Cosine | 24 | 18.03 | 39.80% | 38.58% | 0.2753 | 5.496 / 7.368 |
| Cosine reference | 1,024 | 28.36 | 0.0217% | — | — | 8.068 / 13.585 |
| Gyroid | 8 | 7.32 | 89.98% | 89.42% | 0.7800 | 2.890 / 4.356 |
| Gyroid | 12 | 10.91 | 89.98% | 89.42% | 0.7800 | 4.521 / 10.363 |
| Gyroid | 16 | 14.51 | 88.02% | 87.47% | 0.7607 | 5.894 / 10.485 |
| Gyroid | 24 | 20.92 | 72.12% | 71.57% | 0.6047 | 8.027 / 11.924 |
| Gyroid reference | 1,024 | 58.75 | 0.2267% | — | — | 20.687 / 52.758 |

Cosine starts inside on 58.33% of rays; gyroid starts inside on 50.00%. An exact
zero-set start is emitted at zero depth. The reference is also bounded: grazing
and floating-point progress failures remain unresolved.

Device frame spills, firmware memory, gradient-lighting cost, and verified
supersampling cost require the separate segmented-device experiment. The native
quality result already prevents admitting these low-budget configurations.
