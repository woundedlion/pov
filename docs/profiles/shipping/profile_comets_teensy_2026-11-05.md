# Comets on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Comets`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/comets_ship.log`. Captured 2026-11-04 23:28 America/Los_Angeles on COM4. Replaces profile_comets_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Comets, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 120 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Comets profile 120 16` |

Image size:

```text
teensy_size:   FLASH: code:90376, data:149452, headers:9004   free for files:1782784
teensy_size:    RAM1: variables:313984, code:26968, padding:5800   free for local variables:177536
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `e6d9759dff36661061d1313a1317c92e2828b76fde2c0de30a311c3c379ede9c`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 513–528, root cyc/600 versus measured wall sum differs by **1.135 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Comets                   worst 2 lines/set (0 sets over)  gated 20,376 B  overflow  4,780 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **30.637 ms** (frame 647), **0/1887 spilled** (0.00%). Setup frame 1: 0.860 ms. Complete windows cover frames 1–1888; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 15.368 ms/frame. `cm_draw_trail` averages 14.366 ms/frame over complete windows; its worst window is 20.593 ms/frame (frames 513–528).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 1.023–30.637 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     54.880 ms 32.928 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.3%
    x0.94/frame; 151.205 us/call
  cm_draw_trail                            1.799 ms  1.080 Mcyc   3.3%
    x1.00/frame; 1799.375 us/call
  cm_wipe_rebake                           0.000 ms  0.000 Mcyc   0.0%
    x1.00/frame; 0.096 us/call
  cm_timeline_step                         0.120 ms  0.072 Mcyc   0.2%
    x1.00/frame; 119.839 us/call
  canvas_clear                             0.094 ms  0.056 Mcyc   0.2%
    x1.00/frame; 93.814 us/call
  canvas_buffer_wait                      52.706 ms 31.623 Mcyc  96.0%
    x1.00/frame; 52705.758 us/call
```

Wall min/avg/max: 0.864/54.880/63.490 ms. Window mean render is 2.175 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 641–656)

```text
frame                                     62.065 ms 37.239 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.354 us/call
  cm_draw_trail                           20.043 ms 12.026 Mcyc  32.3%
    x1.00/frame; 20042.866 us/call
  cm_wipe_rebake                           2.026 ms  1.216 Mcyc   3.3%
    x1.00/frame; 2026.277 us/call
  cm_timeline_step                         0.175 ms  0.105 Mcyc   0.3%
    x1.00/frame; 175.266 us/call
  canvas_clear                             0.086 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.799 us/call
  canvas_buffer_wait                      39.570 ms 23.742 Mcyc  63.8%
    x1.00/frame; 39569.740 us/call
```

Wall min/avg/max: 45.660/62.064/77.459 ms. Window mean render is 22.495 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 513–528)

```text
frame                                     62.496 ms 37.497 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 137.786 us/call
  cm_draw_trail                           20.593 ms 12.356 Mcyc  33.0%
    x1.00/frame; 20592.993 us/call
  cm_wipe_rebake                           2.024 ms  1.214 Mcyc   3.2%
    x1.00/frame; 2023.569 us/call
  cm_timeline_step                         0.189 ms  0.114 Mcyc   0.3%
    x1.00/frame; 189.493 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.128 us/call
  canvas_buffer_wait                      39.437 ms 23.662 Mcyc  63.1%
    x1.00/frame; 39437.079 us/call
```

Wall min/avg/max: 58.275/62.495/66.718 ms. Window mean render is 23.058 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 625–640)

```text
frame                                     62.029 ms 37.217 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.547 us/call
  cm_draw_trail                           19.853 ms 11.912 Mcyc  32.0%
    x1.00/frame; 19852.686 us/call
  cm_wipe_rebake                           0.000 ms  0.000 Mcyc   0.0%
    x1.00/frame; 0.094 us/call
  cm_timeline_step                         0.182 ms  0.109 Mcyc   0.3%
    x1.00/frame; 181.683 us/call
  canvas_clear                             0.086 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.751 us/call
  canvas_buffer_wait                      41.729 ms 25.037 Mcyc  67.3%
    x1.00/frame; 41729.054 us/call
```

Wall min/avg/max: 57.097/62.028/66.449 ms. Window mean render is 20.300 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 4  | 10/9 | unavailable | 20.593 | 23.058 | 16.00 | modal-call-count hold |
| 10  | 10/9 | unavailable | 20.080 | 20.511 | 16.12 | modal-call-count hold |
| 5  | 10/9 | unavailable | 20.043 | 22.495 | 16.11 | modal-call-count hold |
| 12  | 8/8 | unavailable | 20.029 | 20.462 | 16.00 | modal-call-count hold |
| 2  | 10/9 | unavailable | 19.602 | 20.030 | 16.01 | modal-call-count hold |
| 8  | 10/9 | unavailable | 19.403 | 19.828 | 15.96 | modal-call-count hold |
| 9  | 10/9 | unavailable | 19.033 | 21.460 | 15.99 | modal-call-count hold |
| 11  | 10/9 | unavailable | 18.069 | 20.610 | 15.98 | modal-call-count hold |
| 6  | 10/9 | unavailable | 15.046 | 17.534 | 15.97 | modal-call-count hold |
| 3  | 10/9 | unavailable | 14.885 | 15.302 | 16.01 | modal-call-count hold |
| 7  | 10/9 | unavailable | 12.511 | 12.914 | 16.01 | modal-call-count hold |
| 1  | 10/9 | unavailable | 12.210 | 12.604 | 16.02 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 5 | 30.637 | 0/160 | 16.00 |
| 11 | 27.578 | 0/160 | 16.00 |
| 10 | 27.121 | 0/160 | 16.00 |
| 4 | 26.260 | 0/160 | 16.00 |
| 9 | 25.898 | 0/160 | 16.00 |
| 2 | 25.074 | 0/160 | 16.00 |
| 12 | 23.461 | 0/129 | 16.00 |
| 8 | 22.355 | 0/160 | 16.00 |
| 6 | 21.126 | 0/160 | 16.00 |
| 3 | 19.184 | 0/160 | 16.00 |
| 7 | 15.800 | 0/160 | 16.00 |
| 1 | 14.833 | 0/158 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1149.56/frame
  min/avg/max 0.563/1.681/21.395 us; cpu 3.10%
  isr_pack         143.69/frame
    min/avg/max 6.176/7.046/10.350 us; cpu 1.62%
  isr_dma_submit   143.69/frame
    min/avg/max 0.620/0.935/1.123 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=997815 us
isr_wake       n=18393 cyc min/avg/max=338/1009/12837 ns min/avg/max=563/1681/21395 total=30946 us cpu=3.10%
isr_pack       n=2299 cyc min/avg/max=3706/4228/6210 ns min/avg/max=6176/7046/10350 total=16202 us cpu=1.62%
isr_dma_submit n=2299 cyc min/avg/max=372/561/674 ns min/avg/max=620/935/1123 total=2151 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.10% leaves approximately 60.562 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `cm_draw_trail` — 20.043 ms/frame, 32.29% of this root window.
2. `cm_wipe_rebake` — 2.026 ms/frame, 3.26% of this root window.
3. `cm_timeline_step` — 0.175 ms/frame, 0.28% of this root window.
4. `pov_preserve_half` — 0.140 ms/frame, 0.23% of this root window.
5. `canvas_clear` — 0.086 ms/frame, 0.14% of this root window.

README cells: peak 🟢 30.637 (12), spilled 🟢 0/1887 (0.00%).

README setup frame 1: 0.860 ms. Captured 2026-11-04 23:28 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Comets  windows=118
  [PASS] captured >=3 windows (118)
  [PASS] every window names one effect and resolution (Comets 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 12, 11 distinct)
  [FAIL] all 12 presets visited (11)
  [PASS] per-frame render is render, not wall (0 of 118 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 118 windows are short)
  [PASS] scope 'cm_draw_trail' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.1 ppm, frames 513-528)

INVALID: Comets
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 12, 11 distinct); [FAIL] all 12 presets visited (11). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Comets`, `HS_PROFILE_WINDOW=16`. `just profile Comets` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
