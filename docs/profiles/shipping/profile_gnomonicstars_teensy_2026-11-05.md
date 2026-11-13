# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GnomonicStars`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/gnomonicstars_ship.log`. Captured 2026-11-04 22:58 America/Los_Angeles on COM4. Replaces profile_gnomonicstars_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GnomonicStars, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh GnomonicStars profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:75992, data:149096, headers:8384   free for files:1798144
teensy_size:    RAM1: variables:314048, code:20840, padding:11928   free for local variables:177472
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `f1de29eae5c02d9f306028d90be117839efc994df872eef3ec0a78d00ff94749`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 193–224, root cyc/600 versus measured wall sum differs by **1.024 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   GnomonicStars            worst 2 lines/set (0 sets over)  gated 20,468 B  overflow  5,928 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **26.237 ms** (frame 206), **0/447 spilled** (0.00%). Setup frame 1: 15.236 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 12.185 ms/frame. `gn_draw_stars` averages 11.856 ms/frame over complete windows; its worst window is 18.237 ms/frame (frames 193–224).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 4.407–26.237 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     59.355 ms 35.613 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x0.97/frame; 142.823 us/call
  gn_draw_stars                           10.167 ms  6.100 Mcyc  17.1%
    gn_star_scan                           8.459 ms  5.075 Mcyc  14.3%
      x251.75/frame; 33.599 us/call
  gn_spiral_build                          0.017 ms  0.010 Mcyc   0.0%
    x0.03/frame; 559.622 us/call
  gn_timeline_step                         0.069 ms  0.041 Mcyc   0.1%
    animation_mobius_step                  0.019 ms  0.011 Mcyc   0.0%
      x1.00/frame; 18.666 us/call
  canvas_clear                             0.094 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.197 us/call
  canvas_buffer_wait                      48.835 ms 29.301 Mcyc  82.3%
    x1.00/frame; 48834.525 us/call
```

Wall min/avg/max: 14.825/59.355/72.410 ms. Window mean render is 10.521 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 193–224)

```text
frame                                     62.494 ms 37.497 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.925 us/call
  gn_draw_stars                           18.237 ms 10.942 Mcyc  29.2%
    gn_star_scan                          16.642 ms  9.985 Mcyc  26.6%
      x238.00/frame; 69.925 us/call
  gn_timeline_step                         0.062 ms  0.037 Mcyc   0.1%
    animation_mobius_step                  0.017 ms  0.010 Mcyc   0.0%
      x1.00/frame; 17.127 us/call
  canvas_clear                             0.089 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.380 us/call
  canvas_buffer_wait                      43.930 ms 26.358 Mcyc  70.3%
    x1.00/frame; 43930.179 us/call
```

Wall min/avg/max: 48.137/62.494/75.807 ms. Window mean render is 18.564 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 97–128)

```text
frame                                     62.371 ms 37.422 Mcyc 100.0%
  pov_preserve_half                        0.146 ms  0.088 Mcyc   0.2%
    x1.00/frame; 146.052 us/call
  gn_draw_stars                            8.670 ms  5.202 Mcyc  13.9%
    gn_star_scan                           6.994 ms  4.197 Mcyc  11.2%
      x228.56/frame; 30.602 us/call
  gn_timeline_step                         0.067 ms  0.040 Mcyc   0.1%
    animation_mobius_step                  0.020 ms  0.012 Mcyc   0.0%
      x1.00/frame; 20.039 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.716 us/call
  canvas_buffer_wait                      53.361 ms 32.017 Mcyc  85.6%
    x1.00/frame; 53361.359 us/call
```

Wall min/avg/max: 56.660/62.370/67.893 ms. Window mean render is 9.009 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.66/frame
  min/avg/max 0.538/1.658/11.835 us; cpu 3.05%
  isr_pack         143.94/frame
    min/avg/max 5.991/6.950/9.751 us; cpu 1.60%
  isr_dma_submit   143.94/frame
    min/avg/max 0.771/0.935/1.056 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999314 us
isr_wake       n=36853 cyc min/avg/max=323/995/7101 ns min/avg/max=538/1658/11835 total=61153 us cpu=3.05%
isr_pack       n=4606 cyc min/avg/max=3595/4170/5851 ns min/avg/max=5991/6950/9751 total=32017 us cpu=1.60%
isr_dma_submit n=4606 cyc min/avg/max=463/561/634 ns min/avg/max=771/935/1056 total=4313 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.05% leaves approximately 60.594 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `gn_draw_stars` — 18.237 ms/frame, 29.18% of this root window.
2. `gn_star_scan` — 16.642 ms/frame, 26.63% of this root window.
3. `pov_preserve_half` — 0.142 ms/frame, 0.23% of this root window.
4. `canvas_clear` — 0.089 ms/frame, 0.14% of this root window.
5. `gn_timeline_step` — 0.062 ms/frame, 0.10% of this root window.

README cells: peak 🟢 26.237, spilled 🟢 0/447 (0.00%).

README setup frame 1: 15.236 ms. Captured 2026-11-04 22:58 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=GnomonicStars  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (GnomonicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'gn_draw_stars' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.0 ppm, frames 193-224)

VALID: GnomonicStars
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`. `just profile GnomonicStars` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
