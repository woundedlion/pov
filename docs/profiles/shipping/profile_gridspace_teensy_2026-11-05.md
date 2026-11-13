# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GridSpace`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/gridspace_ship.log`. Captured 2026-11-04 23:55 America/Los_Angeles on COM4. Replaces profile_gridspace_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridSpace, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh GridSpace profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:100776, data:152412, headers:8956   free for files:1769472
teensy_size:    RAM1: variables:314144, code:28216, padding:4552   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `89d123fbe6eb03e4b4563cf2ec8e40034d60dbd2c8f4fb30e7b365612218ead0`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 897–928, root cyc/600 versus measured wall sum differs by **14.704 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   GridSpace                worst 2 lines/set (0 sets over)  gated 20,402 B  overflow  2,524 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **28.401 ms** (frame 413), **0/927 spilled** (0.00%). Setup frame 1: 41.152 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 23.746 ms/frame. `fx_shader_draw` averages 18.295 ms/frame over complete windows; its worst window is 19.033 ms/frame (frames 897–928).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 21.562–28.401 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     60.461 ms 36.277 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x0.97/frame; 146.617 us/call
  fx_shader_draw                          18.739 ms 11.243 Mcyc  31.0%
    fx_prepare_column                      0.562 ms  0.337 Mcyc   0.9%
      x144.00/frame; 3.903 us/call
  fx_palette_base_bake                     0.185 ms  0.111 Mcyc   0.3%
    fx_palette_rgb_bake                    0.074 ms  0.044 Mcyc   0.1%
      x1.00/frame; 73.858 us/call
  fx_prepare_frame                         3.256 ms  1.953 Mcyc   5.4%
    x1.00/frame; 3255.708 us/call
  fx_advance                               0.979 ms  0.587 Mcyc   1.6%
    palette_cycle_step                     0.890 ms  0.534 Mcyc   1.5%
      x1.00/frame; 890.393 us/call
  fx_timeline_step                         0.109 ms  0.066 Mcyc   0.2%
    x1.00/frame; 109.336 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.805 us/call
  canvas_buffer_wait                      36.723 ms 22.034 Mcyc  60.7%
    x1.00/frame; 36722.850 us/call
```

Wall min/avg/max: 21.967/60.460/65.260 ms. Window mean render is 23.738 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 385–416)

```text
frame                                     62.381 ms 37.429 Mcyc 100.0%
  pov_preserve_half                        0.147 ms  0.088 Mcyc   0.2%
    x1.00/frame; 147.317 us/call
  fx_shader_draw                          18.308 ms 10.985 Mcyc  29.3%
    fx_prepare_column                      0.556 ms  0.334 Mcyc   0.9%
      x144.00/frame; 3.864 us/call
  fx_palette_base_bake                     0.906 ms  0.544 Mcyc   1.5%
    fx_palette_rgb_bake                    0.787 ms  0.472 Mcyc   1.3%
      x1.00/frame; 787.161 us/call
  fx_prepare_frame                         3.334 ms  2.000 Mcyc   5.3%
    x1.00/frame; 3333.778 us/call
  fx_advance                               0.990 ms  0.594 Mcyc   1.6%
    palette_cycle_step                     0.950 ms  0.570 Mcyc   1.5%
      x1.00/frame; 950.235 us/call
  fx_timeline_step                         0.114 ms  0.069 Mcyc   0.2%
    x1.00/frame; 114.391 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.178 us/call
  canvas_buffer_wait                      38.274 ms 22.965 Mcyc  61.4%
    x1.00/frame; 38274.458 us/call
```

Wall min/avg/max: 56.130/62.379/68.644 ms. Window mean render is 24.106 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 897–928)

```text
frame                                     62.376 ms 37.425 Mcyc 100.0%
  pov_preserve_half                        0.148 ms  0.089 Mcyc   0.2%
    x1.00/frame; 148.123 us/call
  fx_shader_draw                          19.033 ms 11.420 Mcyc  30.5%
    fx_prepare_column                      0.554 ms  0.332 Mcyc   0.9%
      x144.00/frame; 3.846 us/call
  fx_palette_base_bake                     1.197 ms  0.718 Mcyc   1.9%
    fx_palette_rgb_bake                    1.082 ms  0.649 Mcyc   1.7%
      x1.00/frame; 1081.860 us/call
  fx_prepare_frame                         3.292 ms  1.975 Mcyc   5.3%
    x1.00/frame; 3292.095 us/call
  fx_advance                               0.996 ms  0.597 Mcyc   1.6%
    palette_cycle_step                     0.958 ms  0.575 Mcyc   1.5%
      x1.00/frame; 958.494 us/call
  fx_timeline_step                         0.119 ms  0.071 Mcyc   0.2%
    x1.00/frame; 118.957 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.084 us/call
  canvas_buffer_wait                      37.279 ms 22.367 Mcyc  59.8%
    x1.00/frame; 37278.807 us/call
```

Wall min/avg/max: 56.425/62.374/68.105 ms. Window mean render is 25.097 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.376 ms 37.426 Mcyc 100.0%
  pov_preserve_half                        0.147 ms  0.088 Mcyc   0.2%
    x1.00/frame; 146.736 us/call
  fx_shader_draw                          18.182 ms 10.909 Mcyc  29.1%
    fx_prepare_column                      0.558 ms  0.335 Mcyc   0.9%
      x144.00/frame; 3.874 us/call
  fx_palette_base_bake                     0.253 ms  0.152 Mcyc   0.4%
    fx_palette_rgb_bake                    0.135 ms  0.081 Mcyc   0.2%
      x1.00/frame; 135.421 us/call
  fx_prepare_frame                         3.272 ms  1.963 Mcyc   5.2%
    x1.00/frame; 3272.457 us/call
  fx_advance                               0.980 ms  0.588 Mcyc   1.6%
    palette_cycle_step                     0.947 ms  0.568 Mcyc   1.5%
      x1.00/frame; 946.514 us/call
  fx_timeline_step                         0.110 ms  0.066 Mcyc   0.2%
    x1.00/frame; 110.236 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.337 us/call
  canvas_buffer_wait                      39.117 ms 23.470 Mcyc  62.7%
    x1.00/frame; 39117.087 us/call
```

Wall min/avg/max: 59.309/62.375/65.485 ms. Window mean render is 23.259 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.91/frame
  min/avg/max 0.540/1.948/24.868 us; cpu 3.59%
  isr_pack         144.00/frame
    min/avg/max 8.206/9.118/11.731 us; cpu 2.10%
  isr_dma_submit   144.00/frame
    min/avg/max 0.665/0.936/1.098 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999737 us
isr_wake       n=36861 cyc min/avg/max=324/1169/14921 ns min/avg/max=540/1948/24868 total=71835 us cpu=3.59%
isr_pack       n=4608 cyc min/avg/max=4924/5471/7039 ns min/avg/max=8206/9118/11731 total=42023 us cpu=2.10%
isr_dma_submit n=4608 cyc min/avg/max=399/562/659 ns min/avg/max=665/936/1098 total=4320 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.59% leaves approximately 60.256 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 18.308 ms/frame, 29.35% of this root window.
2. `fx_prepare_frame` — 3.334 ms/frame, 5.34% of this root window.
3. `fx_advance` — 0.990 ms/frame, 1.59% of this root window.
4. `palette_cycle_step` — 0.950 ms/frame, 1.52% of this root window.
5. `fx_palette_base_bake` — 0.906 ms/frame, 1.45% of this root window.

README cells: peak 🟢 28.401, spilled 🟢 0/927 (0.00%).

README setup frame 1: 41.152 ms. Captured 2026-11-04 23:55 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=GridSpace  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (GridSpace 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (14.7 ppm, frames 897-928)

INVALID: GridSpace
```

Validation flags: [FAIL] root cyc/600 vs wall sum within 5 ppm (14.7 ppm, frames 897-928). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridSpace`, `HS_PROFILE_WINDOW=32`. `just profile GridSpace` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
