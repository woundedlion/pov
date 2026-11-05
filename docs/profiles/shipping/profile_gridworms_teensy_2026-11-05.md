# GridWorms on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GridWorms`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/gridworms_ship.log`. Captured 2026-11-04 23:56 America/Los_Angeles on COM4. Replaces profile_gridworms_teensy_2026-10-31.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridWorms, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh GridWorms profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:94568, data:152572, headers:8860   free for files:1775616
teensy_size:    RAM1: variables:314144, code:28088, padding:4680   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `bbe9a5b2ffbb7b608d5f00f1119e1fd4465dfdba1d40bc2252a245a1ce0b9045`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 97–128, root cyc/600 versus measured wall sum differs by **6.833 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   GridWorms                worst 2 lines/set (0 sets over)  gated 16,386 B  overflow      0 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **37.926 ms** (frame 125), **0/927 spilled** (0.00%). Setup frame 1: 56.868 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 32.488 ms/frame. `fx_shader_draw` averages 29.950 ms/frame over complete windows; its worst window is 33.200 ms/frame (frames 97–128).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 27.101–37.926 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     61.149 ms 36.689 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x0.97/frame; 142.842 us/call
  fx_shader_draw                          29.332 ms 17.599 Mcyc  48.0%
    fx_prepare_column                      1.369 ms  0.821 Mcyc   2.2%
      x144.00/frame; 9.507 us/call
  fx_palette_base_bake                     0.228 ms  0.137 Mcyc   0.4%
    fx_palette_rgb_bake                    0.109 ms  0.065 Mcyc   0.2%
      x1.00/frame; 108.649 us/call
  fx_prepare_frame                         0.409 ms  0.245 Mcyc   0.7%
    x1.00/frame; 408.901 us/call
  fx_advance                               1.028 ms  0.617 Mcyc   1.7%
    palette_cycle_step                     0.931 ms  0.558 Mcyc   1.5%
      x1.00/frame; 930.659 us/call
  fx_timeline_step                         0.115 ms  0.069 Mcyc   0.2%
    x1.00/frame; 115.305 us/call
  canvas_clear                             0.089 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.406 us/call
  canvas_buffer_wait                      29.612 ms 17.767 Mcyc  48.4%
    x1.00/frame; 29612.320 us/call
```

Wall min/avg/max: 29.302/61.148/67.120 ms. Window mean render is 31.536 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 97–128)

```text
frame                                     62.341 ms 37.405 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 142.154 us/call
  fx_shader_draw                          33.200 ms 19.920 Mcyc  53.3%
    fx_prepare_column                      1.445 ms  0.867 Mcyc   2.3%
      x144.00/frame; 10.037 us/call
  fx_palette_base_bake                     0.487 ms  0.292 Mcyc   0.8%
    fx_palette_rgb_bake                    0.368 ms  0.221 Mcyc   0.6%
      x1.00/frame; 367.879 us/call
  fx_prepare_frame                         0.419 ms  0.251 Mcyc   0.7%
    x1.00/frame; 418.907 us/call
  fx_advance                               1.042 ms  0.625 Mcyc   1.7%
    palette_cycle_step                     0.996 ms  0.598 Mcyc   1.6%
      x1.00/frame; 996.428 us/call
  fx_timeline_step                         0.125 ms  0.075 Mcyc   0.2%
    x1.00/frame; 125.197 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.139 us/call
  canvas_buffer_wait                      26.653 ms 15.992 Mcyc  42.8%
    x1.00/frame; 26652.872 us/call
```

Wall min/avg/max: 58.062/62.340/66.460 ms. Window mean render is 35.688 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 641–672)

```text
frame                                     62.300 ms 37.380 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.286 us/call
  fx_shader_draw                          28.302 ms 16.981 Mcyc  45.4%
    fx_prepare_column                      1.451 ms  0.871 Mcyc   2.3%
      x144.00/frame; 10.077 us/call
  fx_palette_base_bake                     0.265 ms  0.159 Mcyc   0.4%
    fx_palette_rgb_bake                    0.145 ms  0.087 Mcyc   0.2%
      x1.00/frame; 145.175 us/call
  fx_prepare_frame                         0.316 ms  0.189 Mcyc   0.5%
    x1.00/frame; 315.798 us/call
  fx_advance                               0.979 ms  0.587 Mcyc   1.6%
    palette_cycle_step                     0.935 ms  0.561 Mcyc   1.5%
      x1.00/frame; 935.087 us/call
  fx_timeline_step                         0.123 ms  0.074 Mcyc   0.2%
    x1.00/frame; 122.736 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 84.829 us/call
  canvas_buffer_wait                      31.896 ms 19.137 Mcyc  51.2%
    x1.00/frame; 31895.639 us/call
```

Wall min/avg/max: 56.332/62.299/68.334 ms. Window mean render is 30.404 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.53/frame
  min/avg/max 0.596/2.050/26.436 us; cpu 3.78%
  isr_pack         144.06/frame
    min/avg/max 6.926/8.766/11.528 us; cpu 2.02%
  isr_dma_submit   144.06/frame
    min/avg/max 0.791/0.956/3.930 us; cpu 0.22%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000822 us
isr_wake       n=36881 cyc min/avg/max=358/1230/15862 ns min/avg/max=596/2050/26436 total=75649 us cpu=3.78%
isr_pack       n=4610 cyc min/avg/max=4156/5260/6917 ns min/avg/max=6926/8766/11528 total=40417 us cpu=2.02%
isr_dma_submit n=4610 cyc min/avg/max=475/574/2358 ns min/avg/max=791/956/3930 total=4412 us cpu=0.22%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.78% leaves approximately 60.137 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 33.200 ms/frame, 53.26% of this root window.
2. `fx_prepare_column` — 1.445 ms/frame, 2.32% of this root window.
3. `fx_advance` — 1.042 ms/frame, 1.67% of this root window.
4. `palette_cycle_step` — 0.996 ms/frame, 1.60% of this root window.
5. `fx_palette_base_bake` — 0.487 ms/frame, 0.78% of this root window.

README cells: peak 🟢 37.926, spilled 🟢 0/927 (0.00%).

README setup frame 1: 56.868 ms. Captured 2026-11-04 23:56 America/Los_Angeles.

Also affected by the 1a3244da8..dc2954193 placement/render fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=GridWorms  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (GridWorms 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (6.8 ppm, frames 97-128)

INVALID: GridWorms
```

Validation flags: [FAIL] root cyc/600 vs wall sum within 5 ppm (6.8 ppm, frames 97-128). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridWorms`, `HS_PROFILE_WINDOW=32`. `just profile GridWorms` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
