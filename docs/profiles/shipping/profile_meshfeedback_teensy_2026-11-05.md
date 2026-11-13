# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/meshfeedback_ship.log`. Captured 2026-11-04 23:32 America/Los_Angeles on COM4. Replaces profile_meshfeedback_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MeshFeedback, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 181 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh MeshFeedback profile 181 16` |

Image size:

```text
teensy_size:   FLASH: code:184144, data:190632, headers:8200   free for files:1648640
teensy_size:    RAM1: variables:314144, code:46344, padding:19192   free for local variables:144608
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `0400b47eb534152c5f0f436077147ea33f0680a5adbdb6e8952635e20a54636f`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 513–528, root cyc/600 versus measured wall sum differs by **2.182 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   MeshFeedback             worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 53,236 B  full-set worst 5
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **53.505 ms** (frame 2651), **0/2863 spilled** (0.00%). Setup frame 1: 6.712 ms. Complete windows cover frames 1–2864; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 44.341 ms/frame. `mf_feedback_flush` averages 36.563 ms/frame over complete windows; its worst window is 40.272 ms/frame (frames 513–528).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 34.031–53.505 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     57.709 ms 34.625 Mcyc 100.0%
  mf_mesh_draw                             6.742 ms  4.045 Mcyc  11.7%
    x1.00/frame; 6742.292 us/call
  mf_feedback_flush                       30.954 ms 18.572 Mcyc  53.6%
    feedback_composite                    23.468 ms 14.081 Mcyc  40.7%
      x0.94/frame; 25032.242 us/call
    feedback_populate                      7.254 ms  4.353 Mcyc  12.6%
      x0.94/frame; 7737.936 us/call
    feedback_litscan                       0.075 ms  0.045 Mcyc   0.1%
      x1.00/frame; 75.275 us/call
  mf_timeline_step                         0.115 ms  0.069 Mcyc   0.2%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.178 us/call
  mf_apply_params                          0.003 ms  0.002 Mcyc   0.0%
    x1.00/frame; 3.143 us/call
  canvas_clear                             0.395 ms  0.237 Mcyc   0.7%
    x1.00/frame; 395.160 us/call
  canvas_buffer_wait                      19.447 ms 11.668 Mcyc  33.7%
    x1.00/frame; 19446.596 us/call
```

Wall min/avg/max: 6.715/57.709/65.181 ms. Window mean render is 38.263 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 2641–2656)

```text
frame                                     62.597 ms 37.558 Mcyc 100.0%
  mf_mesh_draw                             9.419 ms  5.651 Mcyc  15.0%
    x1.00/frame; 9418.695 us/call
  mf_feedback_flush                       32.907 ms 19.744 Mcyc  52.6%
    feedback_composite                    32.264 ms 19.359 Mcyc  51.5%
      x1.00/frame; 32264.226 us/call
    feedback_populate                      0.483 ms  0.290 Mcyc   0.8%
      x1.00/frame; 483.157 us/call
    feedback_litscan                       0.019 ms  0.012 Mcyc   0.0%
      x1.00/frame; 19.254 us/call
  mf_timeline_step                         0.108 ms  0.065 Mcyc   0.2%
    x1.00/frame; 108.090 us/call
  mf_apply_params                          0.043 ms  0.026 Mcyc   0.1%
    x1.00/frame; 43.328 us/call
  canvas_clear                             0.393 ms  0.236 Mcyc   0.6%
    x1.00/frame; 393.257 us/call
  canvas_buffer_wait                      19.678 ms 11.807 Mcyc  31.4%
    x1.00/frame; 19678.365 us/call
```

Wall min/avg/max: 53.968/62.597/74.581 ms. Window mean render is 42.919 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 513–528)

```text
frame                                     62.398 ms 37.439 Mcyc 100.0%
  mf_mesh_draw                             7.730 ms  4.638 Mcyc  12.4%
    x1.00/frame; 7729.998 us/call
  mf_feedback_flush                       40.272 ms 24.163 Mcyc  64.5%
    feedback_composite                    32.320 ms 19.392 Mcyc  51.8%
      x1.00/frame; 32320.462 us/call
    feedback_populate                      7.796 ms  4.678 Mcyc  12.5%
      x1.00/frame; 7795.867 us/call
    feedback_litscan                       0.020 ms  0.012 Mcyc   0.0%
      x1.00/frame; 20.090 us/call
  mf_timeline_step                         0.112 ms  0.067 Mcyc   0.2%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.177 us/call
  mf_apply_params                          0.002 ms  0.001 Mcyc   0.0%
    x1.00/frame; 1.957 us/call
  canvas_clear                             0.398 ms  0.239 Mcyc   0.6%
    x1.00/frame; 398.251 us/call
  canvas_buffer_wait                      13.838 ms  8.303 Mcyc  22.2%
    x1.00/frame; 13837.670 us/call
```

Wall min/avg/max: 62.028/62.398/62.908 ms. Window mean render is 48.561 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Held-preset regime (frames 2129–2144)

```text
frame                                     62.398 ms 37.439 Mcyc 100.0%
  mf_mesh_draw                             8.574 ms  5.144 Mcyc  13.7%
    x1.00/frame; 8573.563 us/call
  mf_feedback_flush                       40.253 ms 24.152 Mcyc  64.5%
    feedback_composite                    32.335 ms 19.401 Mcyc  51.8%
      x1.00/frame; 32334.574 us/call
    feedback_populate                      7.753 ms  4.652 Mcyc  12.4%
      x1.00/frame; 7753.063 us/call
    feedback_litscan                       0.017 ms  0.010 Mcyc   0.0%
      x1.00/frame; 16.580 us/call
  mf_timeline_step                         0.112 ms  0.067 Mcyc   0.2%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.172 us/call
  mf_apply_params                          0.004 ms  0.002 Mcyc   0.0%
    x1.00/frame; 3.624 us/call
  canvas_clear                             0.400 ms  0.240 Mcyc   0.6%
    x1.00/frame; 399.680 us/call
  canvas_buffer_wait                      13.002 ms  7.801 Mcyc  20.8%
    x1.00/frame; 13001.540 us/call
```

Wall min/avg/max: 62.080/62.397/62.696 ms. Window mean render is 49.396 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 1201–1216)

```text
frame                                     62.354 ms 37.412 Mcyc 100.0%
  mf_mesh_draw                             7.924 ms  4.754 Mcyc  12.7%
    x1.00/frame; 7923.549 us/call
  mf_feedback_flush                       39.667 ms 23.800 Mcyc  63.6%
    feedback_composite                    31.751 ms 19.051 Mcyc  50.9%
      x1.00/frame; 31751.189 us/call
    feedback_populate                      7.751 ms  4.651 Mcyc  12.4%
      x1.00/frame; 7751.155 us/call
    feedback_litscan                       0.019 ms  0.012 Mcyc   0.0%
      x1.00/frame; 19.282 us/call
  mf_timeline_step                         0.117 ms  0.070 Mcyc   0.2%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.177 us/call
  mf_apply_params                          0.020 ms  0.012 Mcyc   0.0%
    x1.00/frame; 20.012 us/call
  canvas_clear                             0.398 ms  0.239 Mcyc   0.6%
    x1.00/frame; 397.504 us/call
  canvas_buffer_wait                      14.177 ms  8.506 Mcyc  22.7%
    x1.00/frame; 14176.830 us/call
```

Wall min/avg/max: 60.570/62.353/65.582 ms. Window mean render is 48.177 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 3  | 15/14 | unavailable | 40.272 | 48.561 | 16.03 | modal-call-count hold |
| 9  | 15/14 | unavailable | 40.253 | 49.396 | 16.03 | modal-call-count hold |
| 5  | 15/14 | unavailable | 40.167 | 48.303 | 16.04 | modal-call-count hold |
| 4  | 15/14 | unavailable | 39.515 | 42.847 | 16.04 | modal-call-count hold |
| 6  | 15/14 | unavailable | 39.146 | 47.960 | 16.01 | modal-call-count hold |
| 2  | 16/15 | unavailable | 39.124 | 44.200 | 16.05 | modal-call-count hold |
| 8  | 15/14 | unavailable | 38.754 | 47.951 | 16.02 | modal-call-count hold |
| 7  | 15/14 | unavailable | 36.271 | 40.565 | 16.06 | modal-call-count hold |
| 1  | 15/15 | unavailable | 35.148 | 42.552 | 16.02 | modal-call-count hold |
| 10  | 15/14 | unavailable | 34.171 | 41.047 | 16.05 | modal-call-count hold |
| 12  | 13/13 | unavailable | 33.489 | 45.340 | 16.02 | modal-call-count hold |
| 11  | 15/14 | unavailable | 33.084 | 43.262 | 16.03 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 12 | 53.505 | 0/214 | 16.00 |
| 6 | 50.794 | 0/241 | 16.00 |
| 11 | 50.683 | 0/241 | 16.00 |
| 3 | 49.924 | 0/241 | 16.00 |
| 9 | 49.848 | 0/241 | 16.00 |
| 8 | 49.600 | 0/241 | 16.00 |
| 5 | 48.988 | 0/241 | 16.00 |
| 10 | 48.803 | 0/241 | 16.00 |
| 7 | 45.652 | 0/241 | 16.00 |
| 2 | 44.906 | 0/241 | 16.00 |
| 1 | 43.953 | 0/239 | 16.00 |
| 4 | 43.575 | 0/241 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1153.81/frame
  min/avg/max 0.506/2.038/17.288 us; cpu 3.75%
  isr_pack         144.19/frame
    min/avg/max 7.863/9.631/12.210 us; cpu 2.21%
  isr_dma_submit   144.19/frame
    min/avg/max 0.685/0.938/1.528 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1001478 us
isr_wake       n=18461 cyc min/avg/max=304/1223/10373 ns min/avg/max=506/2038/17288 total=37635 us cpu=3.75%
isr_pack       n=2307 cyc min/avg/max=4718/5779/7326 ns min/avg/max=7863/9631/12210 total=22223 us cpu=2.21%
isr_dma_submit n=2307 cyc min/avg/max=411/563/917 ns min/avg/max=685/938/1528 total=2166 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.75% leaves approximately 60.156 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `mf_feedback_flush` — 32.907 ms/frame, 52.57% of this root window.
2. `feedback_composite` — 32.264 ms/frame, 51.54% of this root window.
3. `mf_mesh_draw` — 9.419 ms/frame, 15.05% of this root window.
4. `feedback_populate` — 0.483 ms/frame, 0.77% of this root window.
5. `canvas_clear` — 0.393 ms/frame, 0.63% of this root window.

README cells: peak 🟢 53.505 (12), spilled 🟢 0/2863 (0.00%).

README setup frame 1: 6.712 ms. Captured 2026-11-04 23:32 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=MeshFeedback  windows=179
  [PASS] captured >=3 windows (179)
  [PASS] every window names one effect and resolution (MeshFeedback 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 12, 11 distinct)
  [FAIL] all 12 presets visited (11)
  [PASS] per-frame render is render, not wall (0 of 179 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 179 windows are short)
  [PASS] scope 'mf_feedback_flush' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (2.2 ppm, frames 513-528)

INVALID: MeshFeedback
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 12, 11 distinct); [FAIL] all 12 presets visited (11). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`. `just profile MeshFeedback` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
