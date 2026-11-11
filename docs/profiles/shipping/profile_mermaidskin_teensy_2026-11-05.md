# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MermaidSkin`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/mermaidskin_ship.log`. Captured 2026-11-05 00:07 America/Los_Angeles on COM4. Replaces profile_mermaidskin_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MermaidSkin, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh MermaidSkin profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:102648, data:151732, headers:8788   free for files:1768448
teensy_size:    RAM1: variables:314144, code:29960, padding:2808   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `9f3304c77c359ce2bd70a28e853a3c6bd207ba50b973793a06995e80e799dfe1`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **0.403 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   MermaidSkin              worst 2 lines/set (0 sets over)  gated 20,478 B  overflow  3,096 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **41.164 ms** (frame 255), **0/447 spilled** (0.00%). Setup frame 1: 64.523 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 36.322 ms/frame. `fx_shader_draw` averages 30.326 ms/frame over complete windows; its worst window is 30.905 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 33.547–41.164 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     61.584 ms 36.951 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.97/frame; 139.644 us/call
  fx_shader_draw                          30.905 ms 18.543 Mcyc  50.2%
    fx_prepare_column                      0.635 ms  0.381 Mcyc   1.0%
      x144.00/frame; 4.407 us/call
  fx_palette_base_bake                     0.234 ms  0.141 Mcyc   0.4%
    fx_palette_rgb_bake                    0.116 ms  0.070 Mcyc   0.2%
      x1.00/frame; 116.211 us/call
  fx_prepare_frame                         3.293 ms  1.976 Mcyc   5.3%
    x1.00/frame; 3292.657 us/call
  fx_advance                               1.208 ms  0.725 Mcyc   2.0%
    palette_cycle_step                     1.109 ms  0.666 Mcyc   1.8%
      x1.00/frame; 1109.476 us/call
  fx_timeline_step                         0.091 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.207 us/call
  canvas_clear                             0.099 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.705 us/call
  canvas_buffer_wait                      25.401 ms 15.241 Mcyc  41.2%
    x1.00/frame; 25401.422 us/call
```

Wall min/avg/max: 33.645/61.584/66.393 ms. Window mean render is 36.183 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 225–256)

```text
frame                                     62.409 ms 37.445 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.963 us/call
  fx_shader_draw                          30.355 ms 18.213 Mcyc  48.6%
    fx_prepare_column                      0.638 ms  0.383 Mcyc   1.0%
      x144.00/frame; 4.434 us/call
  fx_palette_base_bake                     1.222 ms  0.733 Mcyc   2.0%
    fx_palette_rgb_bake                    1.102 ms  0.661 Mcyc   1.8%
      x1.00/frame; 1101.869 us/call
  fx_prepare_frame                         3.366 ms  2.020 Mcyc   5.4%
    x1.00/frame; 3366.382 us/call
  fx_advance                               1.319 ms  0.791 Mcyc   2.1%
    palette_cycle_step                     1.292 ms  0.775 Mcyc   2.1%
      x1.00/frame; 1292.380 us/call
  fx_timeline_step                         0.084 ms  0.051 Mcyc   0.1%
    x1.00/frame; 84.423 us/call
  canvas_clear                             0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 97.097 us/call
  canvas_buffer_wait                      25.626 ms 15.376 Mcyc  41.1%
    x1.00/frame; 25625.916 us/call
```

Wall min/avg/max: 55.713/62.408/69.145 ms. Window mean render is 36.783 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.399 ms 37.439 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.647 us/call
  fx_shader_draw                          30.090 ms 18.054 Mcyc  48.2%
    fx_prepare_column                      0.637 ms  0.382 Mcyc   1.0%
      x144.00/frame; 4.424 us/call
  fx_palette_base_bake                     0.341 ms  0.205 Mcyc   0.5%
    fx_palette_rgb_bake                    0.224 ms  0.134 Mcyc   0.4%
      x1.00/frame; 223.575 us/call
  fx_prepare_frame                         3.303 ms  1.982 Mcyc   5.3%
    x1.00/frame; 3302.875 us/call
  fx_advance                               1.203 ms  0.722 Mcyc   1.9%
    palette_cycle_step                     1.175 ms  0.705 Mcyc   1.9%
      x1.00/frame; 1174.928 us/call
  fx_timeline_step                         0.087 ms  0.052 Mcyc   0.1%
    x1.00/frame; 87.255 us/call
  canvas_clear                             0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 97.221 us/call
  canvas_buffer_wait                      26.934 ms 16.160 Mcyc  43.2%
    x1.00/frame; 26933.551 us/call
```

Wall min/avg/max: 58.567/62.399/66.201 ms. Window mean render is 35.466 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.06/frame
  min/avg/max 0.498/2.001/24.190 us; cpu 3.69%
  isr_pack         144.00/frame
    min/avg/max 8.850/9.540/19.388 us; cpu 2.19%
  isr_dma_submit   144.00/frame
    min/avg/max 0.631/0.933/1.065 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999976 us
isr_wake       n=36866 cyc min/avg/max=299/1201/14514 ns min/avg/max=498/2001/24190 total=73802 us cpu=3.69%
isr_pack       n=4608 cyc min/avg/max=5310/5724/11633 ns min/avg/max=8850/9540/19388 total=43960 us cpu=2.19%
isr_dma_submit n=4608 cyc min/avg/max=379/560/639 ns min/avg/max=631/933/1065 total=4305 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.69% leaves approximately 60.194 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 30.355 ms/frame, 48.64% of this root window.
2. `fx_prepare_frame` — 3.366 ms/frame, 5.39% of this root window.
3. `fx_advance` — 1.319 ms/frame, 2.11% of this root window.
4. `palette_cycle_step` — 1.292 ms/frame, 2.07% of this root window.
5. `fx_palette_base_bake` — 1.222 ms/frame, 1.96% of this root window.

README cells: peak 🟢 41.164, spilled 🟢 0/447 (0.00%).

README setup frame 1: 64.523 ms. Captured 2026-11-05 00:07 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=MermaidSkin  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (MermaidSkin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.4 ppm, frames 1-32)

VALID: MermaidSkin
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MermaidSkin`, `HS_PROFILE_WINDOW=32`. `just profile MermaidSkin` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
