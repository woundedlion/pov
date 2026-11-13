# LiquidWindows on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile LiquidWindows`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/liquidwindows_ship.log`. Captured 2026-11-05 00:17 America/Los_Angeles on COM4. Replaces profile_liquidwindows_teensy_2026-11-04.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | LiquidWindows, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh LiquidWindows profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:95192, data:153144, headers:8688   free for files:1774592
teensy_size:    RAM1: variables:314144, code:24808, padding:7960   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `43f14cdd4d1982d2168101d6a2de1643567f329a580d435e36270048e592d758`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **0.046 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   LiquidWindows            worst 2 lines/set (0 sets over)  gated 20,474 B  overflow  3,488 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **36.585 ms** (frame 257), **0/927 spilled** (0.00%). Setup frame 1: 57.965 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 31.711 ms/frame. `fx_shader_draw` averages 29.226 ms/frame over complete windows; its worst window is 29.916 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 29.487–36.585 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     61.251 ms 36.750 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x0.97/frame; 142.644 us/call
  fx_shader_draw                          29.916 ms 17.949 Mcyc  48.8%
    fx_prepare_column                      1.498 ms  0.899 Mcyc   2.4%
      x144.00/frame; 10.404 us/call
  fx_palette_base_bake                     0.229 ms  0.137 Mcyc   0.4%
    fx_palette_rgb_bake                    0.109 ms  0.065 Mcyc   0.2%
      x1.00/frame; 108.871 us/call
  fx_prepare_frame                         0.401 ms  0.241 Mcyc   0.7%
    x1.00/frame; 401.025 us/call
  fx_advance                               1.016 ms  0.609 Mcyc   1.7%
    palette_cycle_step                     0.940 ms  0.564 Mcyc   1.5%
      x1.00/frame; 940.362 us/call
  fx_timeline_step                         0.018 ms  0.011 Mcyc   0.0%
    x1.00/frame; 17.530 us/call
  canvas_clear                             0.100 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.157 us/call
  canvas_buffer_wait                      29.221 ms 17.533 Mcyc  47.7%
    x1.00/frame; 29221.280 us/call
```

Wall min/avg/max: 29.875/61.250/65.700 ms. Window mean render is 32.029 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 257–288)

```text
frame                                     62.402 ms 37.441 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.235 us/call
  fx_shader_draw                          29.263 ms 17.558 Mcyc  46.9%
    fx_prepare_column                      1.488 ms  0.893 Mcyc   2.4%
      x144.00/frame; 10.333 us/call
  fx_palette_base_bake                     1.086 ms  0.652 Mcyc   1.7%
    fx_palette_rgb_bake                    0.964 ms  0.579 Mcyc   1.5%
      x1.00/frame; 964.306 us/call
  fx_prepare_frame                         0.465 ms  0.279 Mcyc   0.7%
    x1.00/frame; 464.573 us/call
  fx_advance                               0.997 ms  0.598 Mcyc   1.6%
    palette_cycle_step                     0.977 ms  0.586 Mcyc   1.6%
      x1.00/frame; 977.195 us/call
  fx_timeline_step                         0.017 ms  0.010 Mcyc   0.0%
    x1.00/frame; 17.370 us/call
  canvas_clear                             0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 96.890 us/call
  canvas_buffer_wait                      30.136 ms 18.082 Mcyc  48.3%
    x1.00/frame; 30135.975 us/call
```

Wall min/avg/max: 56.675/62.402/68.028 ms. Window mean render is 32.266 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 577–608)

```text
frame                                     62.485 ms 37.491 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.261 us/call
  fx_shader_draw                          29.083 ms 17.450 Mcyc  46.5%
    fx_prepare_column                      1.476 ms  0.885 Mcyc   2.4%
      x144.00/frame; 10.247 us/call
  fx_palette_base_bake                     0.173 ms  0.104 Mcyc   0.3%
    fx_palette_rgb_bake                    0.056 ms  0.034 Mcyc   0.1%
      x1.00/frame; 56.121 us/call
  fx_prepare_frame                         0.311 ms  0.187 Mcyc   0.5%
    x1.00/frame; 311.288 us/call
  fx_advance                               1.005 ms  0.603 Mcyc   1.6%
    palette_cycle_step                     0.985 ms  0.591 Mcyc   1.6%
      x1.00/frame; 984.981 us/call
  fx_timeline_step                         0.018 ms  0.011 Mcyc   0.0%
    x1.00/frame; 17.519 us/call
  canvas_clear                             0.096 ms  0.058 Mcyc   0.2%
    x1.00/frame; 95.983 us/call
  canvas_buffer_wait                      31.439 ms 18.864 Mcyc  50.3%
    x1.00/frame; 31439.327 us/call
```

Wall min/avg/max: 59.605/62.485/65.293 ms. Window mean render is 31.046 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.41/frame
  min/avg/max 0.618/1.995/24.313 us; cpu 3.67%
  isr_pack         144.03/frame
    min/avg/max 8.630/9.381/11.896 us; cpu 2.16%
  isr_dma_submit   144.03/frame
    min/avg/max 0.726/0.938/5.408 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000529 us
isr_wake       n=36877 cyc min/avg/max=371/1197/14588 ns min/avg/max=618/1995/24313 total=73615 us cpu=3.67%
isr_pack       n=4609 cyc min/avg/max=5178/5629/7138 ns min/avg/max=8630/9381/11896 total=43245 us cpu=2.16%
isr_dma_submit n=4609 cyc min/avg/max=436/563/3245 ns min/avg/max=726/938/5408 total=4327 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.67% leaves approximately 60.206 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 29.263 ms/frame, 46.89% of this root window.
2. `fx_prepare_column` — 1.488 ms/frame, 2.38% of this root window.
3. `fx_palette_base_bake` — 1.086 ms/frame, 1.74% of this root window.
4. `fx_advance` — 0.997 ms/frame, 1.60% of this root window.
5. `palette_cycle_step` — 0.977 ms/frame, 1.57% of this root window.

README cells: peak 🟢 36.585, spilled 🟢 0/927 (0.00%).

README setup frame 1: 57.965 ms. Captured 2026-11-05 00:17 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=LiquidWindows  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (LiquidWindows 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.0 ppm, frames 1-32)

VALID: LiquidWindows
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=LiquidWindows`, `HS_PROFILE_WINDOW=32`. `just profile LiquidWindows` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
