# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-10-10, **selective -O3**)

Point-in-time paired comparison of the integer-coordinate direct-AA path.
This supplements [the 2026-10-06 sweep](profile_shapeshifter_teensy_2026-10-06.md); the historical ranking remains a snapshot of that sweep.
Raw captures and attested ELF artifacts are retained in the host scratch directory:

```text
C:/work/temp/review31-fix/aa-baseline.log
C:/work/temp/review31-fix/aa-candidate.log
C:/work/temp/review31-fix/aa-baseline.provenance
C:/work/temp/review31-fix/aa-candidate.provenance
```

## Setup

| | |
|---|---|
| Hardware | Same Teensy 4.0 on COM3 at 600 MHz; flywheel and DMA ISRs live |
| Image | `profile` env, `-Os` plus shipping selective `HS_O3` regions, including direct-AA and scan primitives |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist |
| Method | Baseline then candidate; 155 seconds each, window 16, epoch 1600 revolutions, full nine-preset cycle |
| Reproduce | Set `HS_PROFILE_TREE` to each source checkout, `HS_TEENSY_PORT=COM3`, distinct `HS_PROFILE_OUT`; run `bash tools/profile_one.sh ShapeShifter profile 155 16 -D HS_PROFILE_EPOCH_REVS=1600` |

Captured source commits:

- Baseline: `df2b955eb97618ebd49d97b1c38cd1459b41b40a`.
- Candidate: `eeda8b0a30f5c6bd1336cf5aea69c8f36af931b3`.

The captures ran from clean committed trees. Later parser/report edits do not change these captured sources. Compiler: Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 (20251203).

```text
Baseline profile ELF SHA256
4bad14cae7aa18fa556036641a0e5f41fa7793bdd6bb03feb7f986cf1bfb37d1
Candidate profile ELF SHA256
cea95f9dda8a6ca034f3253c6268624b535f51fc998eedf57a52579e3e41d6c8
Baseline phantasm ELF SHA256
b92b7c09873270c3231c8544f21bef115224167e314a7ffd75a7aeb7b685a704
Candidate phantasm ELF SHA256
46fcf49dc0ba4d01e2589d239881ca33e99d7e5ce99656bb8dce4f7b20a1bcda
```

Candidate profile image: `FLASH: code:89844, data:156856, headers:8272` / `RAM1: variables:315040, code:42776, padding:22760, free:143712` / `RAM2: variables:520064, free:4224`.
The production Phantasm builds move RAM1 code from 173864 to 174040 bytes (+176), with variables fixed at 314784 and FLASH data fixed at 852880 bytes. Neither consumes additional static RAM.

Exactness cross-check: frames 2257–2272 root cycles divided by 600 match measured wall microseconds within **1.8 ppm baseline / 0.5 ppm candidate**. Both parser validations pass all 153 windows, nine distinct presets, wrap to the initial preset, complete per-frame rows, and no epoch reset.

## Frame cadence

Setup frame 1 is excluded. Each pass has 2447 matched frame numbers. Mean render is 19.730 ms baseline / 19.767 ms candidate (+0.19%); peak render is **54.825 / 54.821 ms**, both at frame 2303. Both have **0/2447 (0.00%)** spilled frames.

One display window is 62.5 ms at 480 RPM; one quadrant is 72×144 = 10368 pixels. Every preset holds 16 fps. The `canvas_buffer_wait` scope waits for the next display flip. This one-pass comparison shows unchanged whole-effect cadence and no measurable whole-effect speedup. The integer Circle pole-cap path is a small part of ShapeShifter's mixed workload.

## Phase-by-phase readout

Nine presets run in order, wrap back to entry 1, and finish on entry 2. Values below are candidate per-frame window averages, percentages of the frame root.

### Costliest clean held window (frames 2257–2272)

```text
frame                   62.03 ms  37.22 Mcyc 100.0%
  pov_preserve_half        0.14 ms   0.08 Mcyc   0.2%
  ss_draw_all             43.99 ms  26.39 Mcyc  70.9%
    ss_plot_dispatch        43.65 ms  26.19 Mcyc  70.4%
  ss_timeline_step         0.03 ms   0.02 Mcyc   0.1%
  ss_buffer_wait          17.86 ms  10.71 Mcyc  28.8%
    canvas_clear             0.08 ms   0.05 Mcyc   0.1%
    canvas_buffer_wait      17.77 ms  10.66 Mcyc  28.7%
```

Wall min/avg/max: 49.887 / 62.025 / 73.322 ms. Baseline `ss_draw_all` is 43.992 ms/f; candidate is 43.991 ms/f. The dominant scan/dispatch workload is unchanged at this scale.

### Peak-frame window (frames 2289–2304)

```text
frame                   62.29 ms  37.37 Mcyc 100.0%
  pov_preserve_half        0.14 ms   0.08 Mcyc   0.2%
  ss_draw_all             42.89 ms  25.73 Mcyc  68.9%
    ss_plot_dispatch        42.55 ms  25.53 Mcyc  68.3%
  ss_timeline_step         0.04 ms   0.02 Mcyc   0.1%
  ss_buffer_wait          19.22 ms  11.53 Mcyc  30.8%
    canvas_clear             0.08 ms   0.05 Mcyc   0.1%
    canvas_buffer_wait      19.13 ms  11.48 Mcyc  30.7%
```

Wall min/avg/max: 40.544 / 62.288 / 83.171 ms. Individual wall times include alignment to the flip boundary; peak render remains 54.821 ms, inside the 62.5 ms budget.

### Per-preset table

Costliest modal-call-count clean hold of each preset, from `parse_profile.py ... presets --scope ss_draw_all`; advance-straddling windows excluded. Counts are clean/all windows. All nine presets and the wrap are verified in both captures.

| Entry | Baseline draw ms/f | Candidate draw ms/f | Delta ms/f | Clean/all windows | fps |
|---|--:|--:|--:|--:|--:|
| 1 | 43.992 | 43.991 | -0.001 | 28/30 | 16 |
| 9 | 24.334 | 24.347 | +0.013 | 14/15 | 16 |
| 4 | 23.333 | 23.340 | +0.007 | 14/15 | 16 |
| 8 | 23.210 | 23.334 | +0.124 | 14/14 | 16 |
| 7 | 19.522 | 19.617 | +0.095 | 14/15 | 16 |
| 6 | 16.915 | 16.995 | +0.079 | 14/15 | 16 |
| 2 | 12.757 | 12.823 | +0.066 | 18/19 | 16 |
| 5 | 6.040 | 6.047 | +0.007 | 14/15 | 16 |
| 3 | 5.085 | 5.086 | +0.000 | 14/15 | 16 |

### Per-pixel figures

These timing images do not enable plot-count telemetry or `filter_blend`, so no cycles/blend figure is available. Native framebuffer parity checks cover every integer column in the supported wrapped range, rows -1 through H, seven clips/margins, and alpha values including zero and saturation. The integer counter has its own emitted field; legacy logs parse it as zero.

## Column-ISR / DMA marshaling cost

Candidate peak-frame window, raw readouts (nanoseconds are per call):

```text
isr_wake        1149.2/frame  ns/call 500/1676/11036  cpu 3.09%
isr_pack        143.6/frame  ns/call 6231/6916/9261  cpu 1.59%
isr_dma_submit  143.6/frame  ns/call 643/941/1040  cpu 0.21%
```

`isr_wake` is inclusive of packing and submission; its CPU share alone is the measured flywheel cost. Pack and submit are nested breakdowns. DMA-completion and other interrupts are unmeasured. Wire transfer is asynchronous, rather than an additional synchronous render cost.

## Summary ranking

1. `ss_draw_all`: 43.991 ms/f, about 70.9% of the costliest held frame; `ss_plot_dispatch` accounts for 43.655 ms/f.
2. `pov_preserve_half`: 0.138 ms/f in that window.
3. `ss_timeline_step`: 0.033 ms/f.

The integer overload directly blends one visible pixel and no longer counts ordinary integer samples as clip-boundary splats. Edge-fold rows retain the original splat. Device A/B supports unchanged whole-effect performance; native parity verifies the output. This report makes no native/WASM timing comparison.

README cells: peak 🟢 54.821 (9), spilled 🟢 0/2447 (0.00%).

## Caveats

- All scopes include ISR time because CYCCNT free-runs.
- One pass per image cannot resolve run-to-run uncertainty; the +0.19% aggregate change is not attributed to the optimization.
- Per-pixel scope overhead would affect a count/stall image; neither capture enables those counters. `filter_blend` parenting can hide its subtree when the first parent's call count is zero; no such counter is present here.
- Both runs use shipping selective-O3 direct-AA and scan regions, with no dwell compression. The longer epoch only prevents reinitialization within the full cycle.
- This captures a mixed ShapeShifter workload rather than an isolated integer-plot benchmark. It does not establish a Galaxies performance figure.
- Captured source hashes identify pre-report commits retained in the local branch history; later documentation and parser changes are outside the measured images.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=1600`; use the supported `tools/profile_one.sh` wrapper, which owns the shared device lock for build, flash, and capture.
