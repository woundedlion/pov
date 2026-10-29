# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-29, **selective -O3**)

Peak render over the complete natural shipping lifetime is **57.914 ms on COM3 and 57.823 ms on COM4**, with **0/1,280 spilled frames on each board**. Both measured peaks meet the requested below-58-ms target, by 86 and 177 µs respectively; this is a small measured margin.

The pre-review revision `84d6d90aa6f2b1fbd64d28bdf62b25eec0a804b7` peaks at 78.272/78.136 ms on the same boards. Measured candidate source is `e926bfd4878bc13e9f7205d6715837fe9fe06a89`, based on `c50e6d2dc`. This report replaces the canonical October 23 snapshot; that older source `49bb199f40a655c23983e73ef6739d47ceb7908c` reported 59.437 ms and predates subsequent renderer work.

Raw logs, source/ELF/compiler attestations and retained images are under `C:/work/temp/review-20261029/`: `19-original84-natural-COM{3,4}.log` and `19-final-lattice-natural-COM{3,4}.log`, their `.provenance` files and `artifacts/` directories. Both candidate captures reuse the identical immutable profile and shipping ELFs.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Two Teensy 4.0 boards, COM3/COM4, 600 MHz; flywheel and DMA ISRs live |
| Image | Shipping selective-O3 `profile`, `HS_PROFILE_TARGET=MindSplatter`, window 16; no count, PC, deep or replay flags |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, 288×144 canvas |
| Capture | Supported `tools/profile_one.sh MindSplatter profile 90 16 '-D HS_PROFILE_EPOCH_REVS=640'`, explicit `HS_PROFILE_TREE`, pinned port |
| Lifetime | Normal 640 revolutions at 480 RPM: 80 seconds, including startup and the actual handoff |
| Schedule | All eight authored presets, normal 112-frame holds/48-frame transitions and initial extra hold; deterministic live effect history |
| Change | Shared trail chunks; effect-local ITCM trail/point/lattice entries; cached time-invariant steady shades |

Shipping image: FLASH code 1,009,584 B, data 861,084 B, headers 8,372 B; RAM1 variables 312,736 B, code 173,928 B, padding 22,680 B, free 14,944 B; RAM2 variables 520,064 B, free 4,224 B. The region/layout gate passes. RAM1 code plus padding remains 196,608 B, matching the pre-review bank allocation; variable storage and free space are unchanged.

Profile image: FLASH code 164,152 B, data 550,696 B, headers 9,120 B; RAM1 variables 313,152 B, code 60,600 B, padding 4,936 B, free 145,600 B; RAM2 variables 520,064 B, free 4,224 B.

The extracted first-epoch captures both pass `parse_profile.py ... validate`: 80 complete counter windows, all eight presets and wrap, complete per-frame rows, and root cycles/600 MHz versus wall sum within 0.2 ppm (COM3) and 1.6 ppm (COM4), frames 993–1008.

## Frame cadence

The acceptance metric is the maximum raw frame **render** time, including setup frame 1, through the observed first effect handoff. Render excludes `canvas_buffer_wait`; wall time includes that display synchronization wait. The display budget is 62.5 ms. The final natural peak is frame 1003 on both boards, with 4.586/4.677 ms margin to that display budget.

| Board | Baseline frames | Baseline peak ms | Candidate frames | Candidate peak ms | Peak change | Spills baseline→candidate | Margin below 58 ms |
|---|---:|---:|---:|---:|---:|---:|---:|
| COM3 | 1196 | 78.272 | 1280 | 57.914 | -26.009% | 42→0 | 86 µs |
| COM4 | 1196 | 78.136 | 1280 | 57.823 | -25.997% | 42→0 | 177 µs |

Every raw first-epoch frame is included, including the baseline's trailing 12 rows beyond its last complete counter window. Both baseline runs end after frame 1196 because rendering overruns consume part of their natural 80-second lifetime. The faster candidate reaches frame 1280, including the authored return to preset 1, before the same deadline. The following epoch's rows confirm the actual handoff and are excluded. No artificial lifetime extension or arbitrary 1248/1280-frame cutoff defines the peak.

The first 1196 frames form a matched-history comparison on each board. Every steady frame is faster in this prefix. Setup frame 1 rises from 0.378 to 0.416 ms on COM3 and 0.377 to 0.415 ms on COM4, an increase of 38 µs on each. This remains far below either budget. README cells exclude setup, use the conservative COM3 peak, and retain the separate setup value.

README cells: peak 🟢 57.914 (8), spilled 🟢 0/1279 (0.00%).

## Phase-by-phase readout

### Held preset 7: frames 993–1008, COM3

```
frame                    62.316 ms  37.390 Mcyc 100.0%
  pov_preserve_half       0.136 ms   0.081 Mcyc   0.2%
  msp_draw_particles     49.920 ms  29.952 Mcyc  80.1%
  msp_particle_step       2.954 ms   1.772 Mcyc   4.7%
  msp_timeline_step       0.131 ms   0.079 Mcyc   0.2%
  canvas_clear            0.095 ms   0.057 Mcyc   0.2%
  canvas_buffer_wait      9.054 ms   5.432 Mcyc  14.5%
```

Wall min/avg/max is 58.370/62.316/65.964 ms. Raw render peaks at 57.914 ms. Trail drawing dominates; buffer wait is idle until the display flip.

### Preset 6 transition: frames 913–928, COM3

```
frame                    62.627 ms  37.576 Mcyc 100.0%
  pov_preserve_half       0.135 ms   0.081 Mcyc   0.2%
  msp_draw_particles     38.718 ms  23.231 Mcyc  61.8%
  msp_particle_step       7.742 ms   4.645 Mcyc  12.4%
  msp_timeline_step       0.095 ms   0.057 Mcyc   0.2%
  canvas_clear            0.094 ms   0.056 Mcyc   0.1%
  canvas_buffer_wait     15.814 ms   9.489 Mcyc  25.3%
```

Wall min/avg/max is 52.457/62.627/72.737 ms. The window's raw render peak is 52.795 ms. Transition and particle-state work remain included in the natural-run peak.

### Per-preset matched-history peaks

Both revisions cover frames 1–1196. Peaks include each preset's following transition. All eight per-preset peaks improve on both boards; no heavy window regresses.

| Preset | COM3 baseline→candidate render ms | Change | COM4 baseline→candidate render ms | Change |
|---|---:|---:|---:|---:|
| 1 | 33.104→25.476 | -23.04% | 33.181→25.429 | -23.36% |
| 2 | 64.944→49.541 | -23.72% | 64.586→49.483 | -23.38% |
| 3 | 44.574→32.672 | -26.70% | 44.678→32.673 | -26.87% |
| 4 | 59.279→47.905 | -19.19% | 58.878→47.963 | -18.54% |
| 5 | 46.134→35.806 | -22.39% | 46.026→35.920 | -21.96% |
| 6 | 67.628→52.795 | -21.93% | 67.473→52.857 | -21.66% |
| 7 | 78.272→57.914 | -26.01% | 78.136→57.823 | -26.00% |
| 8 | 41.017→32.401 | -21.01% | 40.707→32.276 | -20.71% |

### Complete candidate natural-run presets

All 1280 frames, including setup and the return to preset 1, are represented here. The scope peak is the greatest clean-hold window mean of `msp_draw_particles`; it is distinct from the raw worst-frame render peak.

| Preset | COM3 peak render ms | COM4 peak render ms | COM3 clean scope peak ms | Spilled/frames per board |
|---|---:|---:|---:|---:|
| 1 | 34.655 | 34.662 | 19.847 | 0/167 |
| 2 | 49.541 | 49.483 | 39.715 | 0/159 |
| 3 | 32.672 | 32.673 | 28.891 | 0/159 |
| 4 | 47.905 | 47.963 | 33.046 | 0/159 |
| 5 | 35.806 | 35.920 | 26.070 | 0/159 |
| 6 | 52.795 | 52.857 | 38.821 | 0/159 |
| 7 | 57.914 | 57.823 | 49.920 | 0/159 |
| 8 | 33.806 | 33.828 | 28.067 | 0/159 |

Preset 1's full natural-run peak, 34.655/34.662 ms, is higher than the baseline's 33.104/33.181 ms because only the candidate reaches eight extra wrap-history frames. In the matched prefix, preset 1 improves to 25.476/25.429 ms. Preset 8 similarly has 159 candidate frames versus 83 baseline frames; the matched-prefix table keeps those workloads comparable.

### Per-pixel figures

No complete blended-pixel scope is available for the modern point/column route, so no cycles-per-pixel estimate is inferred. Native frozen temporal fixtures preserve all 192 sampled framebuffer hashes, 12 aggregate hashes and 12 modern geometry-counter rows. This verifies those fixtures, rather than claiming an on-device framebuffer hash comparison.

## Column-ISR / DMA marshaling cost

Representative held-preset peak window, COM3 frames 993–1008. Values are calls per frame, min/avg/max CPU duration, and share of the measured ISR window:

```
isr_wake         1151.4/frame 0.615/1.917/23.603 us 3.53%
isr_pack          143.9/frame 6.437/8.538/11.347 us 1.96%
isr_dma_submit    143.9/frame 0.730/0.947/9.125 us 0.21%
```

Pack and DMA-submit scopes nest inside wake and must not be added to it. Submit measures CPU setup rather than asynchronous wire time. DMA completion and other IRQs are not fully attributed.

## Summary ranking

1. The complete natural-run peak is below 58 ms on both boards, with a small 86/177-µs measured margin.
2. No frame exceeds 62.5 ms in either 1280-frame candidate lifetime; the pre-review runs spill 42/1196 each.
3. Preset 7 remains the worst regime. Every matched per-preset peak improves; means are secondary to this worst-frame result.
4. ITCM placement and shade reuse retain the temporal path and native fixture output. Region/layout budgets and the existing RAM1 bank allocation remain satisfied.

## Harness

Only the supported wrapper builds, attests, locks, uploads and captures. `HS_PROFILE_EPOCH_REVS=640` restores the normal shipping deadline; the default profiling lifetime would extend the effect beyond that deadline and expose a different workload. Each capture lasts 90 seconds solely to observe the handoff, and only its first epoch enters these results.

This is a two-board live-driver comparison, with deterministic startup, normal preset schedule and matching frame histories through frame 1196. A working frozen on-device replay harness was unavailable; no unrelated harness repair was included. Source/ELF attestations bind each capture to its actual image. The wider revision interval also includes other review changes, so the comparison does not attribute every performance difference solely to the trail refactor.

Native checks include optimized/NDEBUG full-tier MindSplatter, Hopf and Fishbowl suites; exact frozen temporal hashes/counters; and a focused untimed/timed shade test that checks exact RGB, nonzero output, reduced callback count and timed-overload precedence. Scratch negative controls disabling the cache or reversing precedence fail the focused test. The canonical Release suites pass; an exploratory ad hoc O2 build exposed oracle failures that did not reproduce with the canonical configuration. No oracle tolerance or replay pin changed.
