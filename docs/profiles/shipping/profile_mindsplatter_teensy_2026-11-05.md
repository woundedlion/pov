# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/mindsplatter_ship.log`. Captured 2026-11-04 23:16 America/Los_Angeles on COM4. Replaces profile_mindsplatter_teensy_2026-10-30.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MindSplatter, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 80 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh MindSplatter profile 80 16` |

Image size:

```text
teensy_size:   FLASH: code:158680, data:546384, headers:8664   free for files:1317888
teensy_size:    RAM1: variables:314176, code:60824, padding:4712   free for local variables:144576
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `ce7f9d2000dfbe66e2bf417964511e6656cfa102c4fa86a6fc7098cfb80c61f3`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 993–1008, root cyc/600 versus measured wall sum differs by **0.187 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   MindSplatter             worst 2 lines/set (0 sets over)  gated 20,476 B  overflow 42,880 B  full-set worst 5
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **57.896 ms** (frame 1003), **0/1247 spilled** (0.00%). Setup frame 1: 0.389 ms. Complete windows cover frames 1–1248; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 27.230 ms/frame. `msp_draw_particles` averages 24.358 ms/frame over complete windows; its worst window is 49.988 ms/frame (frames 993–1008).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 0.517–57.896 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     54.770 ms 32.862 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.3%
    x0.94/frame; 149.779 us/call
  msp_draw_particles                       0.708 ms  0.425 Mcyc   1.3%
    msp_particle_scan                      0.000 ms  0.000 Mcyc   0.0%
      x0.06/frame; 1.102 us/call
  msp_particle_step                        0.127 ms  0.076 Mcyc   0.2%
    x1.00/frame; 126.952 us/call
  msp_timeline_step                        0.049 ms  0.029 Mcyc   0.1%
    x1.00/frame; 48.809 us/call
  canvas_clear                             0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 95.027 us/call
  canvas_buffer_wait                      53.632 ms 32.179 Mcyc  97.9%
    x1.00/frame; 53631.561 us/call
```

Wall min/avg/max: 0.392/54.769/62.680 ms. Window mean render is 1.138 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 993–1008)

```text
frame                                     62.340 ms 37.404 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.082 Mcyc   0.2%
    x1.00/frame; 136.287 us/call
  msp_draw_particles                      49.988 ms 29.993 Mcyc  80.2%
    x1.00/frame; 49988.061 us/call
  msp_particle_step                        2.975 ms  1.785 Mcyc   4.8%
    x1.00/frame; 2975.100 us/call
  msp_timeline_step                        0.089 ms  0.053 Mcyc   0.1%
    x1.00/frame; 89.063 us/call
  canvas_clear                             0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.379 us/call
  canvas_buffer_wait                       9.034 ms  5.421 Mcyc  14.5%
    x1.00/frame; 9034.430 us/call
```

Wall min/avg/max: 58.316/62.339/66.039 ms. Window mean render is 53.305 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 945–960)

```text
frame                                     62.145 ms 37.287 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 137.440 us/call
  msp_draw_particles                      36.653 ms 21.992 Mcyc  59.0%
    x1.00/frame; 36652.710 us/call
  msp_particle_step                        7.784 ms  4.670 Mcyc  12.5%
    x1.00/frame; 7784.105 us/call
  msp_timeline_step                        0.070 ms  0.042 Mcyc   0.1%
    x1.00/frame; 70.424 us/call
  canvas_clear                             0.087 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.833 us/call
  canvas_buffer_wait                      17.378 ms 10.427 Mcyc  28.0%
    x1.00/frame; 17377.720 us/call
```

Wall min/avg/max: 55.382/62.144/68.097 ms. Window mean render is 44.767 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 7  | 10/9 | unavailable | 49.988 | 53.305 | 16.04 | modal-call-count hold |
| 2  | 10/9 | unavailable | 39.240 | 41.879 | 16.30 | modal-call-count hold |
| 6  | 10/9 | unavailable | 38.772 | 46.858 | 16.07 | modal-call-count hold |
| 4  | 10/9 | unavailable | 33.020 | 36.374 | 16.02 | modal-call-count hold |
| 3  | 10/9 | unavailable | 28.747 | 31.150 | 15.99 | modal-call-count hold |
| 8  | 8/8 | unavailable | 28.002 | 30.426 | 15.97 | modal-call-count hold |
| 5  | 10/9 | unavailable | 25.488 | 27.960 | 16.11 | modal-call-count hold |
| 1  | 10/9 | unavailable | 19.789 | 22.076 | 16.04 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 7 | 57.896 | 0/159 | 16.00 |
| 6 | 52.904 | 0/159 | 16.00 |
| 2 | 49.494 | 0/159 | 16.00 |
| 4 | 47.390 | 0/159 | 16.00 |
| 5 | 35.392 | 0/159 | 16.00 |
| 8 | 32.783 | 0/135 | 16.00 |
| 3 | 32.576 | 0/159 | 16.00 |
| 1 | 25.223 | 0/158 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1155.62/frame
  min/avg/max 0.618/1.868/21.595 us; cpu 3.44%
  isr_pack         144.44/frame
    min/avg/max 6.676/8.290/10.933 us; cpu 1.90%
  isr_dma_submit   144.44/frame
    min/avg/max 0.783/0.931/1.026 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1003037 us
isr_wake       n=18490 cyc min/avg/max=371/1121/12957 ns min/avg/max=618/1868/21595 total=34557 us cpu=3.44%
isr_pack       n=2311 cyc min/avg/max=4006/4974/6560 ns min/avg/max=6676/8290/10933 total=19158 us cpu=1.90%
isr_dma_submit n=2311 cyc min/avg/max=470/559/616 ns min/avg/max=783/931/1026 total=2153 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.44% leaves approximately 60.350 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `msp_draw_particles` — 49.988 ms/frame, 80.19% of this root window.
2. `msp_particle_step` — 2.975 ms/frame, 4.77% of this root window.
3. `pov_preserve_half` — 0.136 ms/frame, 0.22% of this root window.
4. `msp_timeline_step` — 0.089 ms/frame, 0.14% of this root window.
5. `canvas_clear` — 0.088 ms/frame, 0.14% of this root window.

README cells: peak 🟢 57.896 (8), spilled 🟢 0/1247 (0.00%).

README setup frame 1: 0.389 ms. Captured 2026-11-04 23:16 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- The 1,248-frame pass ends before the post-wrap preset-2 entry. Its preset-2
  row is a warm-up visit with a partly filled scene, not steady-load coverage.
  Full-cycle validation needs at least 1,500 frames or
  `-D HS_PROFILE_MSP_WARMUP=1280`; a later full-cycle capture observed spills
  on that entry. The green aggregate covers only this recorded pass.

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=MindSplatter  windows=78
  [PASS] captured >=3 windows (78)
  [PASS] every window names one effect and resolution (MindSplatter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 8, 7 distinct)
  [FAIL] all 8 presets visited (7)
  [PASS] per-frame render is render, not wall (0 of 78 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 78 windows are short)
  [PASS] scope 'msp_draw_particles' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.2 ppm, frames 993-1008)

INVALID: MindSplatter
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 8, 7 distinct); [FAIL] all 8 presets visited (7). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`. `just profile MindSplatter` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
