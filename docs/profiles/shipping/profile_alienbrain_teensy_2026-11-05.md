# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienBrain`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/alienbrain_ship.log`. Captured 2026-11-04 23:48 America/Los_Angeles on COM4. Replaces profile_alienbrain_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienBrain, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 240 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh AlienBrain profile 240 16` |

Image size:

```text
teensy_size:   FLASH: code:103440, data:153176, headers:8600   free for files:1766400
teensy_size:    RAM1: variables:314144, code:28664, padding:4104   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `21e105eacd0496ce1a8c0410c74d0c759d079092306371d8e7137364cc7ebd16`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 881–896, root cyc/600 versus measured wall sum differs by **2.203 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   AlienBrain               worst 2 lines/set (0 sets over)  gated 20,474 B  overflow  3,036 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **32.584 ms** (frame 2046), **0/3807 spilled** (0.00%). Setup frame 1: 53.448 ms. Complete windows cover frames 1–3808; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 27.077 ms/frame. `fx_shader_draw` averages 24.488 ms/frame over complete windows; its worst window is 26.060 ms/frame (frames 881–896).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 23.027–32.584 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     59.539 ms 35.723 Mcyc 100.0%
  pov_preserve_half                        0.129 ms  0.078 Mcyc   0.2%
    x0.94/frame; 138.033 us/call
  fx_shader_draw                          25.967 ms 15.580 Mcyc  43.6%
    fx_prepare_column                      0.639 ms  0.383 Mcyc   1.1%
      x144.00/frame; 4.435 us/call
  fx_palette_base_bake                     0.185 ms  0.111 Mcyc   0.3%
    fx_palette_rgb_bake                    0.059 ms  0.035 Mcyc   0.1%
      x1.00/frame; 58.832 us/call
  fx_prepare_frame                         0.513 ms  0.308 Mcyc   0.9%
    x1.00/frame; 513.445 us/call
  fx_advance                               1.056 ms  0.634 Mcyc   1.8%
    palette_cycle_step                     0.910 ms  0.546 Mcyc   1.5%
      x1.00/frame; 910.048 us/call
  fx_timeline_step                         0.093 ms  0.056 Mcyc   0.2%
    x1.00/frame; 93.093 us/call
  canvas_clear                             0.105 ms  0.063 Mcyc   0.2%
    x1.00/frame; 105.215 us/call
  canvas_buffer_wait                      31.252 ms 18.751 Mcyc  52.5%
    x1.00/frame; 31251.900 us/call
```

Wall min/avg/max: 25.424/59.538/65.146 ms. Window mean render is 28.287 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 2033–2048)

```text
frame                                     62.386 ms 37.431 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.326 us/call
  fx_shader_draw                          25.265 ms 15.159 Mcyc  40.5%
    fx_prepare_column                      1.030 ms  0.618 Mcyc   1.7%
      x144.00/frame; 7.153 us/call
  fx_palette_base_bake                     0.993 ms  0.596 Mcyc   1.6%
    fx_palette_rgb_bake                    0.863 ms  0.518 Mcyc   1.4%
      x1.00/frame; 862.512 us/call
  fx_prepare_frame                         0.406 ms  0.243 Mcyc   0.7%
    x1.00/frame; 405.646 us/call
  fx_advance                               1.039 ms  0.623 Mcyc   1.7%
    palette_cycle_step                     1.012 ms  0.607 Mcyc   1.6%
      x1.00/frame; 1012.204 us/call
  fx_timeline_step                         0.149 ms  0.089 Mcyc   0.2%
    x1.00/frame; 148.628 us/call
  canvas_clear                             0.099 ms  0.060 Mcyc   0.2%
    x1.00/frame; 99.176 us/call
  canvas_buffer_wait                      34.109 ms 20.465 Mcyc  54.7%
    x1.00/frame; 34108.982 us/call
```

Wall min/avg/max: 56.189/62.385/68.624 ms. Window mean render is 28.277 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 881–896)

```text
frame                                     62.381 ms 37.429 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 136.978 us/call
  fx_shader_draw                          26.060 ms 15.636 Mcyc  41.8%
    fx_prepare_column                      1.043 ms  0.626 Mcyc   1.7%
      x144.00/frame; 7.244 us/call
  fx_palette_base_bake                     1.312 ms  0.787 Mcyc   2.1%
    fx_palette_rgb_bake                    1.180 ms  0.708 Mcyc   1.9%
      x1.00/frame; 1180.448 us/call
  fx_prepare_frame                         0.365 ms  0.219 Mcyc   0.6%
    x1.00/frame; 365.200 us/call
  fx_advance                               1.107 ms  0.664 Mcyc   1.8%
    palette_cycle_step                     1.082 ms  0.649 Mcyc   1.7%
      x1.00/frame; 1081.978 us/call
  fx_timeline_step                         0.148 ms  0.089 Mcyc   0.2%
    x1.00/frame; 147.584 us/call
  canvas_clear                             0.100 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.009 us/call
  canvas_buffer_wait                      32.977 ms 19.786 Mcyc  52.9%
    x1.00/frame; 32976.615 us/call
```

Wall min/avg/max: 56.132/62.381/68.537 ms. Window mean render is 29.404 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 593–608)

```text
frame                                     62.586 ms 37.552 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 142.530 us/call
  fx_shader_draw                          24.833 ms 14.900 Mcyc  39.7%
    fx_prepare_column                      0.858 ms  0.515 Mcyc   1.4%
      x144.00/frame; 5.958 us/call
  fx_palette_base_bake                     0.196 ms  0.118 Mcyc   0.3%
    fx_palette_rgb_bake                    0.067 ms  0.040 Mcyc   0.1%
      x1.00/frame; 67.182 us/call
  fx_prepare_frame                         0.390 ms  0.234 Mcyc   0.6%
    x1.00/frame; 389.709 us/call
  fx_advance                               1.101 ms  0.661 Mcyc   1.8%
    palette_cycle_step                     1.069 ms  0.641 Mcyc   1.7%
      x1.00/frame; 1068.869 us/call
  fx_timeline_step                         0.120 ms  0.072 Mcyc   0.2%
    x1.00/frame; 120.313 us/call
  canvas_clear                             0.099 ms  0.060 Mcyc   0.2%
    x1.00/frame; 99.364 us/call
  canvas_buffer_wait                      35.505 ms 21.303 Mcyc  56.7%
    x1.00/frame; 35504.882 us/call
```

Wall min/avg/max: 59.294/62.585/65.656 ms. Window mean render is 27.081 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 2  | 67/66 | unavailable | 26.060 | 29.404 | 16.03 | modal-call-count hold |
| 1  | 38/37 | unavailable | 25.967 | 28.287 | 16.80 | modal-call-count hold |
| 3  | 68/67 | unavailable | 25.310 | 28.438 | 16.07 | modal-call-count hold |
| 4  | 65/65 | unavailable | 25.109 | 27.402 | 16.03 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 3 | 32.584 | 0/1079 | 16.00 |
| 2 | 32.583 | 0/1079 | 16.00 |
| 1 | 32.208 | 0/598 | 16.00 |
| 4 | 32.193 | 0/1051 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.94/frame
  min/avg/max 0.493/1.908/23.216 us; cpu 3.51%
  isr_pack         144.00/frame
    min/avg/max 7.563/8.956/11.425 us; cpu 2.06%
  isr_dma_submit   144.00/frame
    min/avg/max 0.621/0.938/1.033 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999858 us
isr_wake       n=18431 cyc min/avg/max=296/1145/13930 ns min/avg/max=493/1908/23216 total=35189 us cpu=3.51%
isr_pack       n=2304 cyc min/avg/max=4538/5374/6855 ns min/avg/max=7563/8956/11425 total=20637 us cpu=2.06%
isr_dma_submit n=2304 cyc min/avg/max=373/563/620 ns min/avg/max=621/938/1033 total=2165 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.51% leaves approximately 60.306 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 25.265 ms/frame, 40.50% of this root window.
2. `fx_advance` — 1.039 ms/frame, 1.66% of this root window.
3. `fx_prepare_column` — 1.030 ms/frame, 1.65% of this root window.
4. `palette_cycle_step` — 1.012 ms/frame, 1.62% of this root window.
5. `fx_palette_base_bake` — 0.993 ms/frame, 1.59% of this root window.

README cells: peak 🟢 32.584 (4), spilled 🟢 0/3807 (0.00%).

README setup frame 1: 53.448 ms. Captured 2026-11-04 23:48 America/Los_Angeles.

Pullback shader: also affected by the 1a3244da8..dc2954193 cached-phase/color fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=AlienBrain  windows=238
  [PASS] captured >=3 windows (238)
  [PASS] every window names one effect and resolution (AlienBrain 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 4, 3 distinct)
  [FAIL] all 4 presets visited (3)
  [PASS] per-frame render is render, not wall (0 of 238 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 238 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (2.2 ppm, frames 881-896)

INVALID: AlienBrain
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 4, 3 distinct); [FAIL] all 4 presets visited (3). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienBrain`, `HS_PROFILE_WINDOW=16`. `just profile AlienBrain` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
