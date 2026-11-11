# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile SphericalHarmonics`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/sphericalharmonics_ship.log`. Captured 2026-11-04 23:24 America/Los_Angeles on COM4. Replaces profile_sphericalharmonics_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | SphericalHarmonics, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 96 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh SphericalHarmonics profile 96 16` |

Image size:

```text
teensy_size:   FLASH: code:43632, data:145328, headers:8672   free for files:1833984
teensy_size:    RAM1: variables:312896, code:11992, padding:20776   free for local variables:178624
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `08d8d679008dc4f695add890a27cca80cb423db707b914ddcc2e8dc755de0a3a`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1105–1120, root cyc/600 versus measured wall sum differs by **0.302 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   SphericalHarmonics       worst 1 lines/set (0 sets over)  gated  8,004 B  overflow      0 B  full-set worst 1
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **12.951 ms** (frame 1118), **0/1503 spilled** (0.00%). Setup frame 1: 21.277 ms. Complete windows cover frames 1–1504; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 12.109 ms/frame. `sh_rasterize` averages 11.834 ms/frame over complete windows; its worst window is 12.613 ms/frame (frames 1105–1120).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 10.890–12.951 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     56.714 ms 34.029 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x0.94/frame; 146.530 us/call
  sh_rasterize                            11.741 ms  7.044 Mcyc  20.7%
    x1.00/frame; 11740.507 us/call
  sh_timeline_step                         0.036 ms  0.022 Mcyc   0.1%
    x1.00/frame; 35.916 us/call
  canvas_clear                             0.094 ms  0.056 Mcyc   0.2%
    x1.00/frame; 93.818 us/call
  canvas_buffer_wait                      44.695 ms 26.817 Mcyc  78.8%
    x1.00/frame; 44695.150 us/call
```

Wall min/avg/max: 11.355/56.714/62.962 ms. Window mean render is 12.019 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 1105–1120)

```text
frame                                     62.459 ms 37.475 Mcyc 100.0%
  pov_preserve_half                        0.144 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.901 us/call
  sh_rasterize                            12.613 ms  7.568 Mcyc  20.2%
    x1.00/frame; 12612.559 us/call
  sh_timeline_step                         0.036 ms  0.022 Mcyc   0.1%
    x1.00/frame; 35.885 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.020 us/call
  canvas_buffer_wait                      49.571 ms 29.743 Mcyc  79.4%
    x1.00/frame; 49571.027 us/call
```

Wall min/avg/max: 62.227/62.458/62.582 ms. Window mean render is 12.888 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 1137–1152)

```text
frame                                     62.456 ms 37.473 Mcyc 100.0%
  pov_preserve_half                        0.144 ms  0.087 Mcyc   0.2%
    x1.00/frame; 144.256 us/call
  sh_rasterize                            12.511 ms  7.506 Mcyc  20.0%
    x1.00/frame; 12510.800 us/call
  sh_timeline_step                         0.039 ms  0.023 Mcyc   0.1%
    x1.00/frame; 38.709 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 84.907 us/call
  canvas_buffer_wait                      49.667 ms 29.800 Mcyc  79.5%
    x1.00/frame; 49666.888 us/call
```

Wall min/avg/max: 62.364/62.455/62.544 ms. Window mean render is 12.789 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 24  | 4/3 | unavailable | 12.613 | 12.888 | 16.01 | modal-call-count hold |
| 17  | 4/3 | unavailable | 12.574 | 12.850 | 16.01 | modal-call-count hold |
| 16  | 4/3 | unavailable | 12.453 | 12.734 | 16.01 | modal-call-count hold |
| 18  | 4/3 | unavailable | 12.413 | 12.699 | 16.01 | modal-call-count hold |
| 23  | 4/3 | unavailable | 12.373 | 12.659 | 16.01 | modal-call-count hold |
| 22  | 4/3 | unavailable | 12.245 | 12.523 | 16.01 | modal-call-count hold |
| 19  | 4/3 | unavailable | 12.179 | 12.456 | 16.01 | modal-call-count hold |
| 15  | 4/3 | unavailable | 12.100 | 12.377 | 16.01 | modal-call-count hold |
| 10  | 4/3 | unavailable | 12.088 | 12.367 | 16.01 | modal-call-count hold |
| 9  | 4/3 | unavailable | 12.019 | 12.299 | 16.01 | modal-call-count hold |
| 1  | 4/3 | unavailable | 11.982 | 12.268 | 16.01 | modal-call-count hold |
| 14  | 4/3 | unavailable | 11.938 | 12.223 | 16.01 | modal-call-count hold |
| 11  | 4/3 | unavailable | 11.894 | 12.175 | 16.01 | modal-call-count hold |
| 20  | 4/3 | unavailable | 11.861 | 12.142 | 16.01 | modal-call-count hold |
| 21  | 4/3 | unavailable | 11.838 | 12.121 | 16.01 | modal-call-count hold |
| 7  | 4/3 | unavailable | 11.741 | 12.019 | n/a | startup window (includes setup frame 1) |
| 8  | 4/3 | unavailable | 11.672 | 11.945 | 16.01 | modal-call-count hold |
| 12  | 4/3 | unavailable | 11.656 | 11.932 | 16.01 | modal-call-count hold |
| 5  | 4/3 | unavailable | 11.602 | 11.881 | 16.01 | modal-call-count hold |
| 13  | 4/3 | unavailable | 11.593 | 11.872 | 16.01 | modal-call-count hold |
| 4  | 4/3 | unavailable | 11.556 | 11.842 | 16.01 | modal-call-count hold |
| 6  | 2/2 | unavailable | 11.283 | 11.558 | 16.01 | modal-call-count hold |
| 3  | 4/3 | unavailable | 11.069 | 11.351 | 16.02 | modal-call-count hold |
| 2  | 4/3 | unavailable | 11.040 | 11.323 | 16.01 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 24 | 12.951 | 0/64 | 16.00 |
| 17 | 12.892 | 0/64 | 16.00 |
| 1 | 12.802 | 0/64 | 16.00 |
| 18 | 12.800 | 0/64 | 16.00 |
| 16 | 12.795 | 0/64 | 16.00 |
| 23 | 12.708 | 0/64 | 16.00 |
| 22 | 12.611 | 0/64 | 16.00 |
| 19 | 12.536 | 0/64 | 16.00 |
| 20 | 12.530 | 0/64 | 16.00 |
| 15 | 12.437 | 0/64 | 16.00 |
| 10 | 12.398 | 0/64 | 16.00 |
| 9 | 12.392 | 0/64 | 16.00 |
| 14 | 12.281 | 0/64 | 16.00 |
| 11 | 12.267 | 0/64 | 16.00 |
| 21 | 12.221 | 0/64 | 16.00 |
| 12 | 12.126 | 0/64 | 16.00 |
| 2 | 12.074 | 0/64 | 16.00 |
| 13 | 12.049 | 0/64 | 16.00 |
| 8 | 12.018 | 0/64 | 16.00 |
| 4 | 11.936 | 0/64 | 16.00 |
| 5 | 11.916 | 0/64 | 16.00 |
| 7 | 11.726 | 0/62 | 16.00 |
| 6 | 11.689 | 0/33 | 16.00 |
| 3 | 11.471 | 0/64 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.88/frame
  min/avg/max 0.553/1.955/14.558 us; cpu 3.60%
  isr_pack         144.00/frame
    min/avg/max 8.708/9.276/11.940 us; cpu 2.13%
  isr_dma_submit   144.00/frame
    min/avg/max 0.773/0.931/1.006 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999845 us
isr_wake       n=18430 cyc min/avg/max=332/1173/8735 ns min/avg/max=553/1955/14558 total=36052 us cpu=3.60%
isr_pack       n=2304 cyc min/avg/max=5225/5566/7164 ns min/avg/max=8708/9276/11940 total=21376 us cpu=2.13%
isr_dma_submit n=2304 cyc min/avg/max=464/559/604 ns min/avg/max=773/931/1006 total=2148 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.60% leaves approximately 60.250 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `sh_rasterize` — 12.613 ms/frame, 20.19% of this root window.
2. `pov_preserve_half` — 0.144 ms/frame, 0.23% of this root window.
3. `canvas_clear` — 0.085 ms/frame, 0.14% of this root window.
4. `sh_timeline_step` — 0.036 ms/frame, 0.06% of this root window.

README cells: peak 🟢 12.951 (24), spilled 🟢 0/1503 (0.00%).

README setup frame 1: 21.277 ms. Captured 2026-11-04 23:24 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=SphericalHarmonics  windows=94
  [PASS] captured >=3 windows (94)
  [PASS] every window names one effect and resolution (SphericalHarmonics 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 24, 24 distinct)
  [PASS] all 24 presets visited (24)
  [PASS] per-frame render is render, not wall (0 of 94 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 94 windows are short)
  [PASS] scope 'sh_rasterize' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.3 ppm, frames 1105-1120)

VALID: SphericalHarmonics
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=SphericalHarmonics`, `HS_PROFILE_WINDOW=16`. `just profile SphericalHarmonics` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
