# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexOil`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/kaleidoscopehexoil_ship.log`. Captured 2026-11-05 00:26 America/Los_Angeles on COM4. Replaces profile_kaleidoscopehexoil_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeHexOil, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 105 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh KaleidoscopeHexOil profile 105 16` |

Image size:

```text
teensy_size:   FLASH: code:96192, data:156656, headers:8272   free for files:1770496
teensy_size:    RAM1: variables:314144, code:25416, padding:7352   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `4f34fc62db412781a1b0bd8b35a5be8e4f3a1c4fcf530d833df1a336e83cd458`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 817–832, root cyc/600 versus measured wall sum differs by **18.217 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   KaleidoscopeHexOil       worst 2 lines/set (0 sets over)  gated 19,630 B  overflow      0 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **42.739 ms** (frame 888), **0/1647 spilled** (0.00%). Setup frame 1: 64.707 ms. Complete windows cover frames 1–1648; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 35.608 ms/frame. `fx_shader_draw` averages 32.959 ms/frame over complete windows; its worst window is 35.540 ms/frame (frames 817–832).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 30.907–42.739 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     60.719 ms 36.431 Mcyc 100.0%
  pov_preserve_half                        0.134 ms  0.081 Mcyc   0.2%
    x0.94/frame; 143.293 us/call
  fx_shader_draw                          33.809 ms 20.286 Mcyc  55.7%
    fx_prepare_column                      0.687 ms  0.412 Mcyc   1.1%
      x144.00/frame; 4.768 us/call
  fx_palette_base_bake                     0.176 ms  0.106 Mcyc   0.3%
    fx_palette_rgb_bake                    0.067 ms  0.040 Mcyc   0.1%
      x1.00/frame; 67.011 us/call
  fx_prepare_frame                         0.331 ms  0.198 Mcyc   0.5%
    x1.00/frame; 330.801 us/call
  fx_advance                               1.133 ms  0.680 Mcyc   1.9%
    palette_cycle_step                     0.986 ms  0.591 Mcyc   1.6%
      x1.00/frame; 985.611 us/call
  fx_timeline_step                         0.094 ms  0.056 Mcyc   0.2%
    x1.00/frame; 93.553 us/call
  canvas_clear                             0.102 ms  0.061 Mcyc   0.2%
    x1.00/frame; 101.884 us/call
  canvas_buffer_wait                      24.675 ms 14.805 Mcyc  40.6%
    x1.00/frame; 24674.931 us/call
```

Wall min/avg/max: 33.452/60.718/64.710 ms. Window mean render is 36.044 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 881–896)

```text
frame                                     62.394 ms 37.436 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.690 us/call
  fx_shader_draw                          35.225 ms 21.135 Mcyc  56.5%
    fx_prepare_column                      1.533 ms  0.920 Mcyc   2.5%
      x144.00/frame; 10.647 us/call
  fx_palette_base_bake                     1.305 ms  0.783 Mcyc   2.1%
    fx_palette_rgb_bake                    1.185 ms  0.711 Mcyc   1.9%
      x1.00/frame; 1184.800 us/call
  fx_prepare_frame                         0.358 ms  0.215 Mcyc   0.6%
    x1.00/frame; 358.485 us/call
  fx_advance                               1.090 ms  0.654 Mcyc   1.7%
    palette_cycle_step                     1.072 ms  0.643 Mcyc   1.7%
      x1.00/frame; 1071.817 us/call
  fx_timeline_step                         0.169 ms  0.101 Mcyc   0.3%
    x1.00/frame; 168.828 us/call
  canvas_clear                             0.091 ms  0.055 Mcyc   0.1%
    x1.00/frame; 90.990 us/call
  canvas_buffer_wait                      23.815 ms 14.289 Mcyc  38.2%
    x1.00/frame; 23815.174 us/call
```

Wall min/avg/max: 54.167/62.392/70.686 ms. Window mean render is 38.579 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 817–832)

```text
frame                                     62.407 ms 37.444 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 137.141 us/call
  fx_shader_draw                          35.540 ms 21.324 Mcyc  56.9%
    fx_prepare_column                      1.542 ms  0.925 Mcyc   2.5%
      x144.00/frame; 10.709 us/call
  fx_palette_base_bake                     1.257 ms  0.754 Mcyc   2.0%
    fx_palette_rgb_bake                    1.140 ms  0.684 Mcyc   1.8%
      x1.00/frame; 1139.974 us/call
  fx_prepare_frame                         0.365 ms  0.219 Mcyc   0.6%
    x1.00/frame; 364.819 us/call
  fx_advance                               1.132 ms  0.679 Mcyc   1.8%
    palette_cycle_step                     1.112 ms  0.667 Mcyc   1.8%
      x1.00/frame; 1112.131 us/call
  fx_timeline_step                         0.168 ms  0.101 Mcyc   0.3%
    x1.00/frame; 167.737 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.852 us/call
  canvas_buffer_wait                      23.518 ms 14.111 Mcyc  37.7%
    x1.00/frame; 23518.420 us/call
```

Wall min/avg/max: 55.623/62.405/69.220 ms. Window mean render is 38.889 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 593–608)

```text
frame                                     62.580 ms 37.548 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 142.538 us/call
  fx_shader_draw                          31.560 ms 18.936 Mcyc  50.4%
    fx_prepare_column                      1.049 ms  0.629 Mcyc   1.7%
      x144.00/frame; 7.284 us/call
  fx_palette_base_bake                     0.153 ms  0.092 Mcyc   0.2%
    fx_palette_rgb_bake                    0.039 ms  0.024 Mcyc   0.1%
      x1.00/frame; 39.344 us/call
  fx_prepare_frame                         0.383 ms  0.230 Mcyc   0.6%
    x1.00/frame; 382.784 us/call
  fx_advance                               1.146 ms  0.688 Mcyc   1.8%
    palette_cycle_step                     1.122 ms  0.673 Mcyc   1.8%
      x1.00/frame; 1122.424 us/call
  fx_timeline_step                         0.131 ms  0.078 Mcyc   0.2%
    x1.00/frame; 130.534 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.519 us/call
  canvas_buffer_wait                      28.753 ms 17.252 Mcyc  45.9%
    x1.00/frame; 28752.740 us/call
```

Wall min/avg/max: 57.566/62.578/67.210 ms. Window mean render is 33.827 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 2  | 65/65 | unavailable | 35.540 | 38.889 | 16.02 | modal-call-count hold |
| 1  | 38/37 | unavailable | 34.481 | 37.631 | 16.03 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 2 | 42.739 | 0/1049 | 16.00 |
| 1 | 41.825 | 0/598 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.94/frame
  min/avg/max 0.503/1.998/18.701 us; cpu 3.68%
  isr_pack         144.00/frame
    min/avg/max 8.748/9.446/11.768 us; cpu 2.17%
  isr_dma_submit   144.00/frame
    min/avg/max 0.666/0.933/1.033 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999861 us
isr_wake       n=18431 cyc min/avg/max=302/1199/11221 ns min/avg/max=503/1998/18701 total=36856 us cpu=3.68%
isr_pack       n=2304 cyc min/avg/max=5249/5668/7061 ns min/avg/max=8748/9446/11768 total=21767 us cpu=2.17%
isr_dma_submit n=2304 cyc min/avg/max=400/560/620 ns min/avg/max=666/933/1033 total=2151 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.68% leaves approximately 60.200 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 35.225 ms/frame, 56.46% of this root window.
2. `fx_prepare_column` — 1.533 ms/frame, 2.46% of this root window.
3. `fx_palette_base_bake` — 1.305 ms/frame, 2.09% of this root window.
4. `fx_palette_rgb_bake` — 1.185 ms/frame, 1.90% of this root window.
5. `fx_advance` — 1.090 ms/frame, 1.75% of this root window.

README cells: peak 🟢 42.739 (2), spilled 🟢 0/1647 (0.00%).

README setup frame 1: 64.707 ms. Captured 2026-11-05 00:26 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=KaleidoscopeHexOil  windows=103
  [PASS] captured >=3 windows (103)
  [PASS] every window names one effect and resolution (KaleidoscopeHexOil 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 103 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 103 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (18.2 ppm, frames 817-832)

INVALID: KaleidoscopeHexOil
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 2, 1 distinct); [FAIL] all 2 presets visited (1); [FAIL] root cyc/600 vs wall sum within 5 ppm (18.2 ppm, frames 817-832). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexOil`, `HS_PROFILE_WINDOW=16`. `just profile KaleidoscopeHexOil` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
