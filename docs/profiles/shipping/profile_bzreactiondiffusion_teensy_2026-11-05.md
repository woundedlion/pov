# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile BZReactionDiffusion`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/bzreactiondiffusion_ship.log`. Captured 2026-11-04 22:48 America/Los_Angeles on COM4. Replaces profile_bzreactiondiffusion_teensy_2026-10-30.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | BZReactionDiffusion, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh BZReactionDiffusion profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:68792, data:156560, headers:9144   free for files:1797120
teensy_size:    RAM1: variables:314048, code:16200, padding:16568   free for local variables:177472
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `a7086b00d62990823161414c63d5bcd52478fcb36a0e8f296461735b43640048`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **0.634 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   BZReactionDiffusion      worst 2 lines/set (0 sets over)  gated 19,790 B  overflow    740 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **44.574 ms** (frame 596), **0/927 spilled** (0.00%). Setup frame 1: 76.467 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 44.090 ms/frame. `bz_render` averages 43.851 ms/frame over complete windows; its worst window is 44.234 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 41.837–44.574 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.311 ms 37.387 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.97/frame; 139.437 us/call
  bz_render                               44.234 ms 26.541 Mcyc  71.0%
    bz_raster                             39.940 ms 23.964 Mcyc  64.1%
      x2.00/frame; 19969.899 us/call
    bz_physics                             4.199 ms  2.519 Mcyc   6.7%
      x2.00/frame; 2099.305 us/call
    rd_lattice_basis                       0.092 ms  0.055 Mcyc   0.1%
      x1.00/frame; 92.208 us/call
  rd_timeline_step                         0.041 ms  0.025 Mcyc   0.1%
    x1.00/frame; 41.349 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.970 us/call
  canvas_buffer_wait                      17.800 ms 10.680 Mcyc  28.6%
    x1.00/frame; 17800.013 us/call
```

Wall min/avg/max: 42.288/62.310/76.469 ms. Window mean render is 44.511 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 577–608)

```text
frame                                     62.440 ms 37.464 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 140.966 us/call
  bz_render                               43.933 ms 26.360 Mcyc  70.4%
    bz_raster                             39.675 ms 23.805 Mcyc  63.5%
      x2.00/frame; 19837.570 us/call
    bz_physics                             4.163 ms  2.498 Mcyc   6.7%
      x2.00/frame; 2081.639 us/call
    rd_lattice_basis                       0.092 ms  0.055 Mcyc   0.1%
      x1.00/frame; 91.918 us/call
  rd_timeline_step                         0.037 ms  0.022 Mcyc   0.1%
    x1.00/frame; 37.453 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.267 us/call
  canvas_buffer_wait                      18.235 ms 10.941 Mcyc  29.2%
    x1.00/frame; 18234.890 us/call
```

Wall min/avg/max: 61.626/62.439/63.155 ms. Window mean render is 44.205 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 897–928)

```text
frame                                     62.434 ms 37.461 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.220 us/call
  bz_render                               43.724 ms 26.235 Mcyc  70.0%
    bz_raster                             39.453 ms 23.672 Mcyc  63.2%
      x2.00/frame; 19726.672 us/call
    bz_physics                             4.176 ms  2.506 Mcyc   6.7%
      x2.00/frame; 2088.239 us/call
    rd_lattice_basis                       0.092 ms  0.055 Mcyc   0.1%
      x1.00/frame; 91.896 us/call
  rd_timeline_step                         0.037 ms  0.022 Mcyc   0.1%
    x1.00/frame; 37.333 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.043 us/call
  canvas_buffer_wait                      18.440 ms 11.064 Mcyc  29.5%
    x1.00/frame; 18440.276 us/call
```

Wall min/avg/max: 61.958/62.434/62.924 ms. Window mean render is 43.994 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.09/frame
  min/avg/max 0.550/1.990/29.970 us; cpu 3.66%
  isr_pack         144.03/frame
    min/avg/max 8.966/9.540/12.433 us; cpu 2.19%
  isr_dma_submit   144.03/frame
    min/avg/max 0.736/0.938/1.026 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000019 us
isr_wake       n=36867 cyc min/avg/max=330/1194/17982 ns min/avg/max=550/1990/29970 total=73384 us cpu=3.66%
isr_pack       n=4609 cyc min/avg/max=5380/5724/7460 ns min/avg/max=8966/9540/12433 total=43972 us cpu=2.19%
isr_dma_submit n=4609 cyc min/avg/max=442/563/616 ns min/avg/max=736/938/1026 total=4325 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.66% leaves approximately 60.213 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `bz_render` — 43.933 ms/frame, 70.36% of this root window.
2. `bz_raster` — 39.675 ms/frame, 63.54% of this root window.
3. `bz_physics` — 4.163 ms/frame, 6.67% of this root window.
4. `pov_preserve_half` — 0.141 ms/frame, 0.23% of this root window.
5. `rd_lattice_basis` — 0.092 ms/frame, 0.15% of this root window.

README cells: peak 🟢 44.574, spilled 🟢 0/927 (0.00%).

README setup frame 1: 76.467 ms. Captured 2026-11-04 22:48 America/Los_Angeles.

This capture measures the generated hot-code layout against the 11-04 matched sweep; see the sweep scorecard for the delta.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=BZReactionDiffusion  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (BZReactionDiffusion 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'bz_render' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.6 ppm, frames 1-32)

VALID: BZReactionDiffusion
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=BZReactionDiffusion`, `HS_PROFILE_WINDOW=32`. `just profile BZReactionDiffusion` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
