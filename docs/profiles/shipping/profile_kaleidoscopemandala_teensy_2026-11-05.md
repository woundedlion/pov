# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeMandala`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/kaleidoscopemandala_ship.log`. Captured 2026-11-04 23:53 America/Los_Angeles on COM4. Replaces profile_kaleidoscopemandala_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeMandala, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 105 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh KaleidoscopeMandala profile 105 32` |

Image size:

```text
teensy_size:   FLASH: code:103864, data:153616, headers:8760   free for files:1765376
teensy_size:    RAM1: variables:314144, code:28568, padding:4200   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `2ef6803619b7ccc44c333fceb6d585cdd92988bc8077f53c24b4ee5e9e1b705c`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 865–896, root cyc/600 versus measured wall sum differs by **0.564 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   KaleidoscopeMandala      worst 2 lines/set (0 sets over)  gated 20,474 B  overflow  3,500 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **42.737 ms** (frame 896), **0/1631 spilled** (0.00%). Setup frame 1: 62.398 ms. Complete windows cover frames 1–1632; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 34.838 ms/frame. `fx_shader_draw` averages 31.220 ms/frame over complete windows; its worst window is 33.000 ms/frame (frames 865–896).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 29.380–42.737 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     61.390 ms 36.834 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x0.97/frame; 141.584 us/call
  fx_shader_draw                          30.124 ms 18.075 Mcyc  49.1%
    fx_prepare_column                      1.066 ms  0.639 Mcyc   1.7%
      x144.00/frame; 7.400 us/call
  fx_palette_base_bake                     0.220 ms  0.132 Mcyc   0.4%
    fx_palette_rgb_bake                    0.101 ms  0.061 Mcyc   0.2%
      x1.00/frame; 100.945 us/call
  fx_prepare_frame                         0.462 ms  0.277 Mcyc   0.8%
    x1.00/frame; 462.325 us/call
  fx_advance                               1.088 ms  0.653 Mcyc   1.8%
    palette_cycle_step                     0.997 ms  0.598 Mcyc   1.6%
      x1.00/frame; 996.856 us/call
  fx_timeline_step                         0.087 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.981 us/call
  canvas_clear                             0.105 ms  0.063 Mcyc   0.2%
    x1.00/frame; 105.482 us/call
  canvas_buffer_wait                      28.931 ms 17.358 Mcyc  47.1%
    x1.00/frame; 28930.737 us/call
```

Wall min/avg/max: 30.618/61.390/66.033 ms. Window mean render is 32.459 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 865–896)

```text
frame                                     62.380 ms 37.428 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.786 us/call
  fx_shader_draw                          33.000 ms 19.800 Mcyc  52.9%
    fx_prepare_column                      1.437 ms  0.862 Mcyc   2.3%
      x144.00/frame; 9.978 us/call
  fx_palette_base_bake                     1.388 ms  0.833 Mcyc   2.2%
    fx_palette_rgb_bake                    1.263 ms  0.758 Mcyc   2.0%
      x1.00/frame; 1263.166 us/call
  fx_prepare_frame                         3.359 ms  2.016 Mcyc   5.4%
    x1.00/frame; 3359.278 us/call
  fx_advance                               1.083 ms  0.650 Mcyc   1.7%
    palette_cycle_step                     1.052 ms  0.631 Mcyc   1.7%
      x1.00/frame; 1051.548 us/call
  fx_timeline_step                         0.166 ms  0.099 Mcyc   0.3%
    x1.00/frame; 165.780 us/call
  canvas_clear                             0.098 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.127 us/call
  canvas_buffer_wait                      22.938 ms 13.763 Mcyc  36.8%
    x1.00/frame; 22938.242 us/call
```

Wall min/avg/max: 55.753/62.379/69.041 ms. Window mean render is 39.441 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 577–608)

```text
frame                                     62.591 ms 37.554 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.497 us/call
  fx_shader_draw                          32.367 ms 19.420 Mcyc  51.7%
    fx_prepare_column                      1.183 ms  0.710 Mcyc   1.9%
      x144.00/frame; 8.214 us/call
  fx_palette_base_bake                     0.260 ms  0.156 Mcyc   0.4%
    fx_palette_rgb_bake                    0.135 ms  0.081 Mcyc   0.2%
      x1.00/frame; 135.119 us/call
  fx_prepare_frame                         1.165 ms  0.699 Mcyc   1.9%
    x1.00/frame; 1164.949 us/call
  fx_advance                               1.089 ms  0.653 Mcyc   1.7%
    palette_cycle_step                     1.057 ms  0.634 Mcyc   1.7%
      x1.00/frame; 1056.727 us/call
  fx_timeline_step                         0.103 ms  0.062 Mcyc   0.2%
    x1.00/frame; 102.780 us/call
  canvas_clear                             0.100 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.064 us/call
  canvas_buffer_wait                      27.142 ms 16.285 Mcyc  43.4%
    x1.00/frame; 27141.876 us/call
```

Wall min/avg/max: 58.022/62.590/67.498 ms. Window mean render is 35.449 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 2  | 32/32 | unavailable | 33.000 | 39.441 | 16.03 | modal-call-count hold |
| 1  | 19/18 | unavailable | 32.437 | 35.070 | 16.02 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 2 | 42.737 | 0/1033 | 16.00 |
| 1 | 39.682 | 0/598 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.34/frame
  min/avg/max 0.535/2.008/31.368 us; cpu 3.70%
  isr_pack         144.03/frame
    min/avg/max 8.725/9.413/12.395 us; cpu 2.16%
  isr_dma_submit   144.03/frame
    min/avg/max 0.610/0.941/2.558 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000487 us
isr_wake       n=36875 cyc min/avg/max=321/1205/18821 ns min/avg/max=535/2008/31368 total=74063 us cpu=3.70%
isr_pack       n=4609 cyc min/avg/max=5235/5648/7437 ns min/avg/max=8725/9413/12395 total=43392 us cpu=2.16%
isr_dma_submit n=4609 cyc min/avg/max=366/565/1535 ns min/avg/max=610/941/2558 total=4347 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.70% leaves approximately 60.188 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 33.000 ms/frame, 52.90% of this root window.
2. `fx_prepare_frame` — 3.359 ms/frame, 5.39% of this root window.
3. `fx_prepare_column` — 1.437 ms/frame, 2.30% of this root window.
4. `fx_palette_base_bake` — 1.388 ms/frame, 2.22% of this root window.
5. `fx_palette_rgb_bake` — 1.263 ms/frame, 2.02% of this root window.

README cells: peak 🟢 42.737 (2), spilled 🟢 0/1631 (0.00%).

README setup frame 1: 62.398 ms. Captured 2026-11-04 23:53 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=KaleidoscopeMandala  windows=51
  [PASS] captured >=3 windows (51)
  [PASS] every window names one effect and resolution (KaleidoscopeMandala 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 51 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 51 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.6 ppm, frames 865-896)

INVALID: KaleidoscopeMandala
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 2, 1 distinct); [FAIL] all 2 presets visited (1). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeMandala`, `HS_PROFILE_WINDOW=32`. `just profile KaleidoscopeMandala` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
