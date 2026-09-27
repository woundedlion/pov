# Periodic surface device experiment (2026-09-26)

Cosine and gyroid remain **experimental and excluded from firmware admission**.
Every 8–24-query configuration has substantial unfinished searches and misses
the 62.5 ms display deadline. A black fallback for unfinished rays is an
experimental visualization, not evidence of empty geometry.

## Method and provenance

The [experiment target](../../tests/profile_spherical_experiment.h) runs on
COM3, Teensy 4.0 at 600 MHz, through the existing profile wrapper and real
four-segment 288×144 driver. Three separate shipping-selective-O3 images measure
depth coloring, analytic-gradient lighting, and gradient lighting with four
verified subrays per pixel. No field approximation or enlarged hit tolerance is
introduced. Query budgets include all verification probes; step and refinement
limits are 1,024, with position tolerance equal to period / 10,000.

Each image cycles cosine then gyroid at 8, 12, 16, and 24 queries. Each case
visits all 36 combinations of periods 0.7, 1.5, and 3; isovalues −0.5, 0, and
0.5; and camera centers (0,0,0), (0.23,0.41,−0.17), (0.5,0.5,0.5), and
(−0.37,0.19,0.73), scaled by the period. Radial offset and near distance are
zero; far distance is four periods. The two-sided first-boundary search can
start in either region.

Each reported row uses the first complete 36-frame live block for its case.
Frame 1 is the setup draw and is excluded. The initial cosine/8 block is
therefore replaced by its next complete block after the cycle wraps. Timing
uses exact per-frame render telemetry, including per-ray diagnostic integer
increments and one approximately 170-byte diagnostic serial line per frame.
ISR work remains included. These are instrumented experimental costs, not a
claim about an optimized final renderer.

The active display quadrant is 144×72. The existing shader includes its
one-pixel clip margin, producing **10,658 actual rays** (146×73) for a live
single-sample frame and 42,632 for four samples. The setup draw covers 288×73.
The driver selects the left or right column half according to its current
display window; spills can change that selection. The native quality sweep
evaluates both margin-expanded column halves for every position. Its quality/image-error measurements
and these device costs have separate sample denominators and are not identical
ray-for-ray captures. Separate lighting captures can also select different
column halves; their difference includes that phase variation.

The frozen source is baseline commit
79a93c6f8f39a0538794803b6e45ee7c7710cdab plus the four explicit source files in
the [SHA-256 manifest](evidence/spherical_periodic_device_2026-09-26/source_manifest.json).
Exact snapshots are retained with the evidence: [contracts](evidence/spherical_periodic_device_2026-09-26/contract.h.txt),
[marcher](evidence/spherical_periodic_device_2026-09-26/march.h.txt),
[periodic definitions](evidence/spherical_periodic_device_2026-09-26/periodic_surface.h.txt),
and [experiment](evidence/spherical_periodic_device_2026-09-26/profile_spherical_experiment.h.txt).
This snapshot precedes the additional footprint validation and bit-based
periodic-parameter finite guards. Every captured parameter is finite and the
footprint is zero; the numerical search and query formulas match the final
experiment contracts. Timings remain attributed to the retained snapshot.
The only other source overlay is one include of the experiment header immediately
after the phantasm target include in the profile sketch. The standard wrapper
records compiler/ABI, flag and ELF hashes, builds and gates the full Phantasm
image, and holds the shared device lock through upload and capture.
Environment dumps are retained as JSON strings with their original bytes and
SHA-256 values; readable build logs have trailing whitespace removed.

## Measurements

| Capture | Raw log | Per-case summary | ELF/compiler provenance | Build and size |
| --- | --- | --- | --- | --- |
| Depth | [log](evidence/spherical_periodic_device_2026-09-26/depth.txt) | [JSON](evidence/spherical_periodic_device_2026-09-26/depth_summary.json) | [provenance](evidence/spherical_periodic_device_2026-09-26/depth.provenance) | [build](evidence/spherical_periodic_device_2026-09-26/depth_build.txt) |
| Gradient | [log](evidence/spherical_periodic_device_2026-09-26/gradient.txt) | [JSON](evidence/spherical_periodic_device_2026-09-26/gradient_summary.json) | [provenance](evidence/spherical_periodic_device_2026-09-26/gradient.provenance) | [build](evidence/spherical_periodic_device_2026-09-26/gradient_build.txt) |
| Gradient, four samples | [log](evidence/spherical_periodic_device_2026-09-26/supersampled.txt) | [JSON](evidence/spherical_periodic_device_2026-09-26/supersampled_summary.json) | [provenance](evidence/spherical_periodic_device_2026-09-26/supersampled.provenance) | [build](evidence/spherical_periodic_device_2026-09-26/supersampled_build.txt) |

All three untouched logs pass parser validation and cover the eight-case wrap.
Root cycles / 600 versus wall microseconds agree within 0.4 ppm in the checked
253–288 window. The captures contain 11, 12, and 13 complete windows,
respectively; each selected row below contains 36 live frames. Single-sample
rows contain 383,688 rays and four-sample rows contain 1,534,752 subrays.

| Surface | Queries | Depth mean / peak ms | Gradient mean / peak ms | Depth unfinished | Depth spills |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cosine | 8 | 124.808 / 128.786 | 126.285 / 129.452 | 90.20% | 36/36 |
| Cosine | 12 | 171.239 / 181.201 | 173.866 / 182.389 | 76.88% | 36/36 |
| Cosine | 16 | 205.353 / 235.257 | 205.490 / 233.774 | 60.94% | 36/36 |
| Cosine | 24 | 265.876 / 321.356 | 273.553 / 325.642 | 41.92% | 36/36 |
| Gyroid | 8 | 159.030 / 184.360 | 160.145 / 184.610 | 88.52% | 33/36 |
| Gyroid | 12 | 231.338 / 266.248 | 233.450 / 266.362 | 89.35% | 33/36 |
| Gyroid | 16 | 303.385 / 347.899 | 304.346 / 348.001 | 89.13% | 33/36 |
| Gyroid | 24 | 441.809 / 509.867 | 443.773 / 510.227 | 82.43% | 33/36 |

| Surface | Queries | Four-sample gradient mean / peak ms | Unfinished subrays | Spills |
| --- | ---: | ---: | ---: | ---: |
| Cosine | 8 | 506.397 / 518.864 | 90.14% | 36/36 |
| Cosine | 12 | 685.006 / 730.758 | 72.23% | 36/36 |
| Cosine | 16 | 845.199 / 935.738 | 61.82% | 36/36 |
| Cosine | 24 | 1104.373 / 1303.819 | 44.10% | 36/36 |
| Gyroid | 8 | 649.657 / 739.356 | 90.18% | 36/36 |
| Gyroid | 12 | 922.525 / 1066.282 | 88.52% | 36/36 |
| Gyroid | 16 | 1222.431 / 1391.728 | 89.31% | 36/36 |
| Gyroid | 24 | 1778.277 / 2041.768 | 82.35% | 36/36 |

Unfinished combines unresolved, exhausted, and invalid terminal statuses.
The raw summaries preserve those categories separately; no invalid query was
observed. All searches terminate at their declared work limits. Missed deadlines
and unfinished searches independently prevent admission.

## Memory

| Image | FLASH code | FLASH data | RAM1 code | RAM1 variables | RAM2 variables |
| --- | ---: | ---: | ---: | ---: | ---: |
| Depth, one sample | 35,152 B | 143,864 B | 17,464 B | 314,912 B | 520,064 B |
| Gradient, one sample | 35,808 B | 143,864 B | 17,464 B | 314,912 B | 520,064 B |
| Gradient, four samples | 36,272 B | 143,864 B | 17,464 B | 314,912 B | 520,064 B |

Each image leaves 176,608 B for RAM1 locals and 4,224 B for RAM2 allocation.
Gradient lighting adds 656 B of FLASH code to this experiment; four-sample
filtering adds another 464 B. RAM1 code, variables, and FLASH data do not grow.

The experiment adds seven 32-bit aggregate counters and a 32-bit frame cursor
to the base effect: 32 bytes of persistent state, with no allocated geometry or
per-ray heap storage. Value-type sizes are Ray 32 bytes, TraceResult 72 bytes,
TraceLimits 24 bytes, and periodic surface 20 bytes. These are object sizes,
not an assertion about optimized stack high-water usage. The profile image
contains both periodic variants and all eight quality cases.

The unchanged baseline full-roster Phantasm gate passed at RAM1 code 195,304 B,
RAM1 variables 314,784 B, FLASH data 725,460 B, and FLASH code 501,392 B.
Its remaining ITCM bank headroom was 1,304 B. These experimental image sizes do
not imply that adding every pattern to the full roster meets that budget.

Fresh fixed-preset HyperLattice baselines at that source on COM3 have 1,087 live
frames each: cubic peak 40.984 ms and hypercube peak 50.174 ms, both with zero
spills. The [baseline summary](evidence/spherical_periodic_device_2026-09-26/baseline_summary.json)
retains those measurements and the 96-passed/one-skipped native baseline
(120 smoke frames, no failures). Setup draws are excluded.

## Reproduction

Use the frozen baseline and source overlay above. Add the experiment header
include to the profile sketch after its existing phantasm target include; this
temporary hook is not part of the shipping target. Run each capture through the
supported wrapper, pinning the intended tree and board:

```sh
HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh SphericalExperiment profile 110 36
HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh SphericalExperiment profile 115 36 "-D HS_SPHERICAL_PROFILE_GRADIENT=1"
HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh SphericalExperiment profile 430 36 "-D HS_SPHERICAL_PROFILE_GRADIENT=1 -D HS_SPHERICAL_PROFILE_SAMPLES=4"
```

Preserve the capture between commands because the wrapper uses the same log
name. Validate every untouched capture with the profile parser and confirm
all eight markers plus the wrap. Join each `spherical counts` line's frame
number to the corresponding per-frame render record. Select complete blocks
of 36 positions, excluding frame 1; sum query/status totals and divide by the
actual sample count. Spills count frames with render time above 62,500 μs.

The native high-budget reference and image differences are documented in the
[companion native experiment](spherical_periodic_native_2026-09-26.md). A torus, warped torus, lattice
volume, and triangular framework have native architectural demonstrations;
they receive no device-admission claim from these periodic measurements.
