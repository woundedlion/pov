# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-10-07, **selective -O3**)

Experimental A–B–A characterization. The bounded-clock candidate is not an applied shipping fix. This supplement compares a clean instrument-only baseline with eight bounded double-precision channel phases; it does not replace the canonical GnomonicStars shipping profile.

The measured `animation_mobius_step` scope costs **13.591323 us/frame** in B1 versus **14.263004 us/frame** pooled across A1/A2: **-0.671682 us (-4.709%)**. This is a net scope measurement including eight trig calls, scope overhead and interrupt time. It does not isolate the added clock arithmetic. The complete render is **+9.056578 us/frame** across post-setup frames; the experiment provides no whole-frame performance gain claim.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Same Teensy 4.0 on COM3, 600 MHz; flywheel and DMA ISRs live |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Image | `profile`, shipping `-Os` + newlib-nano, existing selective `HS_O3` regions |
| Hot path | Terminal filter pipeline and annotated math helpers retain selective -O3; Mobius step has no added -O3 annotation |
| Effect | GnomonicStars 288×144, steady single-entry playlist; 600 stars, four sides, default warp speed 0.035 |
| Method | Supported `profile_one.sh` wrapper, 70 s per capture, 32-frame dumps; no dwell compression |
| Baseline A | `17390a7a1eef599b82b1fc13ca531c9156f69f06`, one instrumentation scope over source `e7f19b0d1ef46966f85165cd48cb71b5b34c0d78` |
| Candidate B | `2619870dd2cc9528433e55e7f2e2b1ff687414ea`, same baseline plus bounded Mobius channels and regression tests |
| Compiler | Arm GNU 15.2.Rel1, GCC 15.2.1 20251203 |

Capture order was A1 → B1 → A2. A2's first upload attempt returned Teensy Loader busy while another board was active; a normal supported-wrapper retry completed. The successful A2 capture alone enters this comparison. Each source tree was clean. The native candidate animation suite passed 25,823 assertions; removing channel wrapping produced 76 intended failures (reported by the capture owner).

Raw captures and their sidecars are under `C:/work/temp/finding4-profile-20261007/`:

| Run | Raw capture | Profile ELF SHA-256 | Phantasm ELF SHA-256 |
|---|---|---|---|
| A1 | `gnomonicstars_baseline1.log` | `74030844a7f676537ad65be6578f5ad575046eef6a77debdc5c13f2928826c5a` | `7b60b22ccd2edc55ea24101c5c96308843838df202bfaa039d0dadf7a1462616` |
| B1 | `gnomonicstars_candidate1.log` | `df1195c27e731f8a00fb9bc8de46509b86637fe198c9c5d5a9c76d35d88fced5` | `d13e5dd9f4d0f0fbfbea2a3f629bfe3f7b16cc56736f006d0489d4352efbb825` |
| A2 | `gnomonicstars_baseline2.log` | `f6e64565a03a067359f6b36ad00854e45db9a55b00ece1f7894b4caf85d98dd0` | `28c97f053397e7f36d2220f22bcf00840fd85a66eb14c9af6313542903fa4f5c` |

The `.provenance` sidecars retain the exact source, compiler, environment-dump hashes and saved ELF artifact paths for each successful capture. Different ELF hashes between repeats are retained rather than treated as identical artifacts.

### Image and allocation cost

| Image/run | FLASH code | FLASH data | FLASH headers | RAM1 code | RAM1 variables | RAM1 padding | RAM1 free | RAM2 variables/free |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| profile A1 | 54980 | 153136 | 8968 | 25352 | 315104 | 7416 | 176416 | 520064/4224 |
| profile B1 | 55148 | 153160 | 8776 | 25496 | 315104 | 7272 | 176416 | 520064/4224 |
| profile A2 | 54980 | 153136 | 8968 | 25352 | 315104 | 7416 | 176416 | 520064/4224 |
| phantasm A1 | 635588 | 846284 | 9068 | 177656 | 314784 | 18952 | 12896 | 520064/4224 |
| phantasm B1 | 635748 | 846316 | 8876 | 177816 | 314784 | 18792 | 12896 | 520064/4224 |
| phantasm A2 | 635588 | 846284 | 9068 | 177656 | 314784 | 18952 | 12896 | 520064/4224 |

B1 adds 144 bytes of profile RAM1 code and 24 bytes of FLASH data; shipping Phantasm adds 160 bytes of RAM1 code and 32 bytes of FLASH data. RAM1 variables, RAM2 usage and RAM1 free space remain equal in these images because the code increments fit existing ITCM padding.

The device animation object remains 100 bytes within the 112-byte timeline slot. Its phase state separately consumes 64 bytes from the persistent arena, plus up to seven alignment bytes per construction. The GnomonicStars persistent-footprint formula includes that allocation. This dynamic cost is not visible as increased static RAM1 variables.

### Exactness cross-check

All three captures pass `parse_profile.py ... validate`: 34 complete windows, one effect/resolution, monotonic frame numbers, no epoch reset, complete per-frame telemetry, and render separated from display-sync wait. The richest window is frames 193–224 in each run:

| Run | Root cycles | Root cycles / 600 (us) | Measured wall sum (us) | Difference (ppm) |
|---|---:|---:|---:|---:|
| A1 | 1201156686 | 2001927.810000 | 2001931 | 1.593462 |
| B1 | 1201166410 | 2001944.016667 | 2001948 | 1.989729 |
| A2 | 1201149725 | 2001916.208333 | 2001920 | 1.894015 |

## Frame cadence

All three passes hold the 62.5 ms render budget with zero spilled frames. These per-frame figures exclude setup frame 1 and include frames 2–1088; clock-scope comparisons below conservatively exclude the entire startup window 1–32.

| Run | Post-setup frames | Mean render (ms) | Exact peak render (ms) | Peak frame | Spilled |
|---|---:|---:|---:|---:|---|
| A1 | 1087 | 10.241865 | 21.841 | 210 | 0/1087 |
| B1 | 1087 | 10.250793 | 21.840 | 210 | 0/1087 |
| A2 | 1087 | 10.241608 | 21.844 | 210 | 0/1087 |

A2−A1 mean render is -0.256670 us. B1−pooled A is +9.056578 us. On the stricter matched clock subset (frames 33–1088), mean render is 10.266815/10.275677/10.266568 ms for A1/B1/A2, a B1−pooled A delta of +8.985322 us. The peak remains near 21.84 ms, with approximately 40.66 ms of budget margin.

Small phase differences can alter transformed geometry, raster work and interrupt placement. Render or star-scan deltas cannot be assigned directly to the clock arithmetic. `canvas_buffer_wait` is intentional display-sync idle.

## Phase-by-phase readout

This is a steady effect, with no preset schedule. Match exact frame-window endpoints across runs, discard frames 1–32, and use all 33 common windows (1056 frames and Mobius calls per run, frames 33–1088). The weighted scope cost is `sum(cycles) / sum(calls) / 600`, rather than an unweighted average of rounded printed microseconds.

| Run | Weighted Mobius (us/call) | Window min / median / max (us/call) | Timeline (us/frame) | Draw stars (ms/frame) |
|---|---:|---|---:|---:|
| A1 | 14.286545 | 10.130260 / 14.626615 / 19.338385 | 43.749053 | 9.991596 |
| B1 | 13.591323 | 12.335833 / 13.560104 / 14.877292 | 42.725008 | 10.001769 |
| A2 | 14.239463 | 9.988646 / 14.465938 / 19.245677 | 43.708400 | 9.991490 |

A2−A1 Mobius cost is -0.047082 us/call. B1 is -0.671682 us/call relative to the pooled baselines. This repeat spread is an observed repeat check, not a confidence interval; only one B run was captured. Window variation includes phase-dependent trig argument reduction and ISR placement.

| Matched frames | Windows | A1 Mobius (us) | B1 Mobius (us) | A2 Mobius (us) | B1−mean A (us) |
|---|---:|---:|---:|---:|---:|
| 33–64 | 1 | 12.150990 | 12.747552 | 12.026302 | +0.658906 |
| 33–192 | 5 | 13.795167 | 13.660323 | 13.756562 | -0.115542 |
| 193–512 | 10 | 13.140599 | 13.509224 | 13.143667 | +0.367091 |
| 513–1088 | 18 | 15.059676 | 13.617766 | 14.982378 | -1.403261 |

The first matched window is slightly more expensive with bounded channels. The net advantage appears in later windows; it is not a fixed improvement on every frame. Of the 33 window comparisons, 11 are positive and 22 negative. The largest positive delta is **+2.338932 us/call at frames 193–224**; the largest negative delta is **-5.071406 us/call at frames 993–1024**.

Raw worst-positive scope rows: A1 = 194501 cycles / 32 calls; B1 = 238049 / 32; A2 = 191782 / 32. Dividing each cycle total by 32 calls and 600 cycles/us gives 10.130260, 12.398385, 9.988646 us/call respectively. The complete matched-window table and raw counter trees are retained in `comparison.json` and `analysis.json`.

### Matched-window counter trees

Each tree below averages frames 33–1088. Indentation mirrors the recorded nesting. Percentages are relative to the immediate parent; idle dominates the root.

A1:

```text
frame                         62.457464 ms 37.4745 Mcyc 100.00%
  pov_preserve_half           0.145320 ms  0.0872 Mcyc   0.23%
  gn_draw_stars               9.991596 ms  5.9950 Mcyc  16.00%
    gn_star_scan              9.279612 ms  5.5678 Mcyc  92.87%
      filter_blend            0.338899 ms  0.2033 Mcyc   3.65%
  gn_timeline_step            0.043749 ms  0.0262 Mcyc   0.07%
    animation_mobius_step     0.014287 ms  0.0086 Mcyc  32.66%
  canvas_clear                0.084403 ms  0.0506 Mcyc   0.14%
  canvas_buffer_wait          52.191302 ms 31.3148 Mcyc  83.56%
```

B1:

```text
frame                         62.457644 ms 37.4746 Mcyc 100.00%
  pov_preserve_half           0.145115 ms  0.0871 Mcyc   0.23%
  gn_draw_stars               10.001769 ms  6.0011 Mcyc  16.01%
    gn_star_scan              9.290290 ms  5.5742 Mcyc  92.89%
      filter_blend            0.339028 ms  0.2034 Mcyc   3.65%
  gn_timeline_step            0.042725 ms  0.0256 Mcyc   0.07%
    animation_mobius_step     0.013591 ms  0.0082 Mcyc  31.81%
  canvas_clear                0.084312 ms  0.0506 Mcyc   0.13%
  canvas_buffer_wait          52.182603 ms 31.3096 Mcyc  83.55%
```

A2:

```text
frame                         62.457464 ms 37.4745 Mcyc 100.00%
  pov_preserve_half           0.145220 ms  0.0871 Mcyc   0.23%
  gn_draw_stars               9.991490 ms  5.9949 Mcyc  16.00%
    gn_star_scan              9.279512 ms  5.5677 Mcyc  92.87%
      filter_blend            0.338378 ms  0.2030 Mcyc   3.65%
  gn_timeline_step            0.043708 ms  0.0262 Mcyc   0.07%
    animation_mobius_step     0.014239 ms  0.0085 Mcyc  32.58%
  canvas_clear                0.084413 ms  0.0506 Mcyc   0.14%
  canvas_buffer_wait          52.191525 ms 31.3149 Mcyc  83.56%
```

### Per-pixel figures

The segmented effect renders one quadrant, approximately 10,368 pixels, with 600 star scans per frame. Blend calls count writes, not unique pixel coverage.

| Run | Blend calls/frame | Cycles/blend | Star-scan cycles/blend |
|---|---:|---:|---:|
| A1 | 3173.201705 | 64.080163 | 1754.621246 |
| B1 | 3173.192235 | 64.104866 | 1756.645498 |
| A2 | 3173.201705 | 63.981666 | 1754.602442 |

## Column-ISR / DMA marshaling cost

Aggregate the dedicated ISR timestamp intervals corresponding to the matched readouts. Those intervals include logging time and differ from the frame-scope intervals. Means below use summed printed ISR totals and counts; window totals have integer-microsecond rounding.

| Run | ISR | Calls/s | Min / mean / max (us) | CPU share |
|---|---|---:|---|---:|
| A1 | `isr_wake` | 18433.183 | 0.650000 / 1.698998 / 17.475000 | 3.131794% |
| A1 | `isr_pack` | 2304.000 | 6.238333 / 6.739927 / 9.681667 | 1.552879% |
| A1 | `isr_dma_submit` | 2304.000 | 0.608333 / 0.933685 / 1.155000 | 0.215121% |
| B1 | `isr_wake` | 18433.184 | 0.541667 / 1.638898 / 17.948333 | 3.021010% |
| B1 | `isr_pack` | 2304.000 | 6.238333 / 6.738690 / 9.518333 | 1.552594% |
| B1 | `isr_dma_submit` | 2304.000 | 0.618333 / 0.938702 / 1.126667 | 0.216277% |
| A2 | `isr_wake` | 18433.185 | 0.653333 / 1.698877 / 20.675000 | 3.131571% |
| A2 | `isr_pack` | 2304.000 | 6.238333 / 6.739105 / 9.576667 | 1.552690% |
| A2 | `isr_dma_submit` | 2304.000 | 0.611667 / 0.933770 / 1.238333 | 0.215141% |

`isr_wake` includes pack and submit; never add their CPU shares to wake. Wake share is slightly lower in B1, while pack and submit are close across runs. No ISR cost is subtracted from measured render or Mobius scopes. Submit is substantially cheaper than pack; asynchronous wire-transfer duration and DMA completion cost are outside these accumulators.

## Summary ranking

1. Star scan/raster work dominates render, approximately 9.28–9.29 ms/frame in matched windows.
2. The complete timeline is approximately 42.73–43.75 us/frame; the measured Mobius scope is approximately 13.59–14.29 us/frame.
3. Bounded channels reduce the net Mobius scope average in this capture, with a small positive early-window cost and larger negative later-window differences. Complete-render performance and cadence remain essentially unchanged at this scale.

GnomonicStars is a legitimate case with the most proposed clock work: eight independent phase advances and eight sinf/cosf evaluations each frame. Noise and NoiseProduct each integrate one live time axis and remain unchanged in this candidate. A future bounded scalar triangle clock for those animations was neither implemented nor measured here; these numbers cannot predict its cost or field semantics.

## Caveats

- CYCCNT free-runs: all profile scopes absorb interrupt time. The direct Mobius scope includes clock arithmetic, eight trig calls and profile overhead. Identical instrumentation placement is common to A and B; actual interrupt placement and execution cost can still differ.
- Bounded trig arguments can reduce library argument-reduction work, plausibly explaining the lower later-window cost. This is an inference, not a separate measured breakdown.
- The candidate preserves the original binary32 channel frequencies, seed-derived offsets, sin/cos choice, scale/base coefficients and live finite speed/scale setters. Negative and zero finite speeds remain accepted; non-finite setter values retain the previous good value.
- Float-max speed times the largest frequency remains finite after promotion to double. The ordinary shipping slider uses speed 0–1, so phase increments need at most one wrap; captures do not measure the large-speed fmod fallback.
- Extremely small negative steps can round a wrapped phase to exactly the period, an equivalent endpoint angle. The prototype guarantees bounded phase magnitude, not mathematically exact advancement for every finite subnormal speed.
- Eight doubles are persistent-arena owned. Copies share the phase block: independently stepping copies advances shared state. The normal pinned GnomonicStars timeline owns one active animation; general copy/respawn use requires attention to shared phase ownership, arena lifetime and allocation until reset.
- Phase integration changes numerical trajectories relative to unbounded binary32 accumulation. Tiny raster differences and ISR placement changes prevent treating the complete render delta as pure clock overhead. A 70-second pass does not reproduce days-long float freeze; native tests exercise large advances and wrapping.
- `filter_blend` is nested beneath `gn_star_scan`; its calls approximate blended writes. Existing selective -O3 annotations remain active on the device build. No cycling, dwell compression or epoch crossing was used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`. The capture owner used the supported wrapper separately with `HS_PROFILE_TREE` pointing to each clean baseline/candidate worktree and `HS_TEENSY_PORT=COM3`:

```sh
bash tools/profile_one.sh GnomonicStars profile 70 32
```

Read-only validation and comparison:

```sh
python tools/parse_profile.py <completed-log> validate
python C:/work/temp/finding4-profile-20261007/analyze.py
python C:/work/temp/finding4-profile-20261007/compare.py
python C:/work/temp/finding4-profile-20261007/build_comparison_report.py
```

Only completed, successfully validated captures were analyzed. The [canonical baseline report](profile_gnomonicstars_teensy_2026-10-07.md) retains the original float clock. This supplemental comparison does not adopt the experimental bounded implementation.
