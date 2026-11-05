# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeFlowers`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/kaleidoscopeflowers_ship.log`. Captured 2026-11-05 00:39 America/Los_Angeles on COM4. Replaces profile_kaleidoscopeflowers_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeFlowers, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 173 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh KaleidoscopeFlowers profile 173 16` |

Image size:

```text
teensy_size:   FLASH: code:103832, data:153336, headers:9072   free for files:1765376
teensy_size:    RAM1: variables:314144, code:28664, padding:4104   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `1d7d521a6fd227a4e78d463e8652487a55de5fc111ddcacb8bc9cfe64ca7a654`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1953–1968, root cyc/600 versus measured wall sum differs by **13.187 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   KaleidoscopeFlowers      worst 2 lines/set (0 sets over)  gated 20,474 B  overflow  3,464 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **36.868 ms** (frame 1966), **0/2735 spilled** (0.00%). Setup frame 1: 57.037 ms. Complete windows cover frames 1–2736; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 30.483 ms/frame. `fx_shader_draw` averages 27.893 ms/frame over complete windows; its worst window is 30.437 ms/frame (frames 1953–1968).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 26.716–36.868 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     59.866 ms 35.920 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.94/frame; 144.193 us/call
  fx_shader_draw                          27.717 ms 16.630 Mcyc  46.3%
    fx_prepare_column                      0.752 ms  0.451 Mcyc   1.3%
      x144.00/frame; 5.219 us/call
  fx_palette_base_bake                     0.183 ms  0.110 Mcyc   0.3%
    fx_palette_rgb_bake                    0.068 ms  0.041 Mcyc   0.1%
      x1.00/frame; 68.155 us/call
  fx_prepare_frame                         0.592 ms  0.355 Mcyc   1.0%
    x1.00/frame; 592.169 us/call
  fx_advance                               1.129 ms  0.677 Mcyc   1.9%
    palette_cycle_step                     0.972 ms  0.583 Mcyc   1.6%
      x1.00/frame; 972.228 us/call
  fx_timeline_step                         0.080 ms  0.048 Mcyc   0.1%
    x1.00/frame; 79.694 us/call
  canvas_clear                             0.108 ms  0.065 Mcyc   0.2%
    x1.00/frame; 108.037 us/call
  canvas_buffer_wait                      29.690 ms 17.814 Mcyc  49.6%
    x1.00/frame; 29689.840 us/call
```

Wall min/avg/max: 27.308/59.865/65.601 ms. Window mean render is 30.176 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 1953–1968)

```text
frame                                     62.595 ms 37.557 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 137.612 us/call
  fx_shader_draw                          30.437 ms 18.262 Mcyc  48.6%
    fx_prepare_column                      1.785 ms  1.071 Mcyc   2.9%
      x144.00/frame; 12.398 us/call
  fx_palette_base_bake                     1.370 ms  0.822 Mcyc   2.2%
    fx_palette_rgb_bake                    1.243 ms  0.746 Mcyc   2.0%
      x1.00/frame; 1242.743 us/call
  fx_prepare_frame                         0.449 ms  0.270 Mcyc   0.7%
    x1.00/frame; 449.413 us/call
  fx_advance                               1.024 ms  0.614 Mcyc   1.6%
    palette_cycle_step                     0.998 ms  0.599 Mcyc   1.6%
      x1.00/frame; 997.816 us/call
  fx_timeline_step                         0.168 ms  0.101 Mcyc   0.3%
    x1.00/frame; 168.061 us/call
  canvas_clear                             0.100 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.174 us/call
  canvas_buffer_wait                      28.756 ms 17.254 Mcyc  45.9%
    x1.00/frame; 28756.248 us/call
```

Wall min/avg/max: 57.570/62.594/67.283 ms. Window mean render is 33.839 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Held-preset regime (frames 1969–1984)

```text
frame                                     62.369 ms 37.422 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.399 us/call
  fx_shader_draw                          30.420 ms 18.252 Mcyc  48.8%
    fx_prepare_column                      1.769 ms  1.061 Mcyc   2.8%
      x144.00/frame; 12.286 us/call
  fx_palette_base_bake                     1.512 ms  0.907 Mcyc   2.4%
    fx_palette_rgb_bake                    1.388 ms  0.833 Mcyc   2.2%
      x1.00/frame; 1387.572 us/call
  fx_prepare_frame                         0.437 ms  0.262 Mcyc   0.7%
    x1.00/frame; 436.673 us/call
  fx_advance                               1.035 ms  0.621 Mcyc   1.7%
    palette_cycle_step                     1.010 ms  0.606 Mcyc   1.6%
      x1.00/frame; 1009.935 us/call
  fx_timeline_step                         0.166 ms  0.099 Mcyc   0.3%
    x1.00/frame; 165.560 us/call
  canvas_clear                             0.099 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.934 us/call
  canvas_buffer_wait                      28.409 ms 17.046 Mcyc  45.5%
    x1.00/frame; 28409.187 us/call
```

Wall min/avg/max: 57.368/62.368/67.203 ms. Window mean render is 33.960 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 1665–1680)

```text
frame                                     62.465 ms 37.479 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.354 us/call
  fx_shader_draw                          28.524 ms 17.114 Mcyc  45.7%
    fx_prepare_column                      0.822 ms  0.493 Mcyc   1.3%
      x144.00/frame; 5.710 us/call
  fx_palette_base_bake                     0.921 ms  0.553 Mcyc   1.5%
    fx_palette_rgb_bake                    0.800 ms  0.480 Mcyc   1.3%
      x1.00/frame; 800.113 us/call
  fx_prepare_frame                         0.459 ms  0.275 Mcyc   0.7%
    x1.00/frame; 458.737 us/call
  fx_advance                               0.982 ms  0.589 Mcyc   1.6%
    palette_cycle_step                     0.950 ms  0.570 Mcyc   1.5%
      x1.00/frame; 950.321 us/call
  fx_timeline_step                         0.082 ms  0.049 Mcyc   0.1%
    x1.00/frame; 82.446 us/call
  canvas_clear                             0.099 ms  0.060 Mcyc   0.2%
    x1.00/frame; 99.220 us/call
  canvas_buffer_wait                      31.045 ms 18.627 Mcyc  49.7%
    x1.00/frame; 31044.769 us/call
```

Wall min/avg/max: 60.196/62.463/64.629 ms. Window mean render is 31.420 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 3  | 66/66 | unavailable | 30.437 | 33.839 | 15.98 | modal-call-count hold |
| 2  | 67/66 | unavailable | 29.442 | 32.136 | 16.02 | modal-call-count hold |
| 1  | 38/37 | unavailable | 29.032 | 31.166 | 16.02 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 3 | 36.868 | 0/1058 | 16.00 |
| 1 | 35.350 | 0/598 | 16.00 |
| 2 | 33.993 | 0/1079 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.56/frame
  min/avg/max 0.546/2.000/17.981 us; cpu 3.68%
  isr_pack         143.94/frame
    min/avg/max 8.858/9.486/12.008 us; cpu 2.18%
  isr_dma_submit   143.94/frame
    min/avg/max 0.626/0.938/3.748 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999592 us
isr_wake       n=18425 cyc min/avg/max=328/1200/10789 ns min/avg/max=546/2000/17981 total=36852 us cpu=3.68%
isr_pack       n=2303 cyc min/avg/max=5315/5692/7205 ns min/avg/max=8858/9486/12008 total=21851 us cpu=2.18%
isr_dma_submit n=2303 cyc min/avg/max=376/563/2249 ns min/avg/max=626/938/3748 total=2163 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.68% leaves approximately 60.200 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 30.437 ms/frame, 48.63% of this root window.
2. `fx_prepare_column` — 1.785 ms/frame, 2.85% of this root window.
3. `fx_palette_base_bake` — 1.370 ms/frame, 2.19% of this root window.
4. `fx_palette_rgb_bake` — 1.243 ms/frame, 1.99% of this root window.
5. `fx_advance` — 1.024 ms/frame, 1.64% of this root window.

README cells: peak 🟢 36.868 (3), spilled 🟢 0/2735 (0.00%).

README setup frame 1: 57.037 ms. Captured 2026-11-05 00:39 America/Los_Angeles.

Pullback shader: also affected by the 1a3244da8..dc2954193 cached-phase/color fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=KaleidoscopeFlowers  windows=171
  [PASS] captured >=3 windows (171)
  [PASS] every window names one effect and resolution (KaleidoscopeFlowers 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 3, 2 distinct)
  [FAIL] all 3 presets visited (2)
  [PASS] per-frame render is render, not wall (0 of 171 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 171 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (13.2 ppm, frames 1953-1968)

INVALID: KaleidoscopeFlowers
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 3, 2 distinct); [FAIL] all 3 presets visited (2); [FAIL] root cyc/600 vs wall sum within 5 ppm (13.2 ppm, frames 1953-1968). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeFlowers`, `HS_PROFILE_WINDOW=16`. `just profile KaleidoscopeFlowers` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
