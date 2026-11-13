# Galaxies on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Galaxies`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/galaxies_ship.log`. Captured 2026-11-04 22:44 America/Los_Angeles on COM4. Replaces profile_galaxies_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Galaxies, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Galaxies profile 60 16` |

Image size:

```text
teensy_size:   FLASH: code:91800, data:149656, headers:8400   free for files:1781760
teensy_size:    RAM1: variables:313984, code:23608, padding:9160   free for local variables:177536
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `8978407eefd1dce060a87aed0ffbac079d4c71335231c1077088b1f218335da2`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 865–880, root cyc/600 versus measured wall sum differs by **0.467 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Galaxies                 worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 20,246 B  full-set worst 3
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **45.126 ms** (frame 883), **0/927 spilled** (0.00%). Setup frame 1: 1.482 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 24.948 ms/frame. `gx_particle_step` averages 15.528 ms/frame over complete windows; its worst window is 27.123 ms/frame (frames 865–880).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 1.254–45.126 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     54.822 ms 32.893 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.3%
    x0.94/frame; 150.501 us/call
  gx_draw_particles                        0.830 ms  0.498 Mcyc   1.5%
    x1.00/frame; 829.871 us/call
  gx_particle_step                         0.427 ms  0.256 Mcyc   0.8%
    x1.00/frame; 427.187 us/call
  gx_timeline_step                         0.077 ms  0.046 Mcyc   0.1%
    x1.00/frame; 77.037 us/call
  canvas_clear                             0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.618 us/call
  canvas_buffer_wait                      53.225 ms 31.935 Mcyc  97.1%
    x1.00/frame; 53224.654 us/call
```

Wall min/avg/max: 1.348/54.822/62.690 ms. Window mean render is 1.597 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 881–896)

```text
frame                                     62.279 ms 37.367 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 145.344 us/call
  gx_draw_particles                       16.454 ms  9.872 Mcyc  26.4%
    x1.00/frame; 16453.646 us/call
  gx_particle_step                        27.101 ms 16.261 Mcyc  43.5%
    x1.00/frame; 27101.289 us/call
  gx_timeline_step                         0.087 ms  0.052 Mcyc   0.1%
    x1.00/frame; 87.425 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.293 us/call
  canvas_buffer_wait                      18.377 ms 11.026 Mcyc  29.5%
    x1.00/frame; 18377.475 us/call
```

Wall min/avg/max: 60.762/62.278/63.783 ms. Window mean render is 43.901 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 865–880)

```text
frame                                     62.389 ms 37.433 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 145.279 us/call
  gx_draw_particles                       16.879 ms 10.128 Mcyc  27.1%
    x1.00/frame; 16879.314 us/call
  gx_particle_step                        27.123 ms 16.274 Mcyc  43.5%
    x1.00/frame; 27123.358 us/call
  gx_timeline_step                         0.089 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.535 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.481 us/call
  canvas_buffer_wait                      18.038 ms 10.823 Mcyc  28.9%
    x1.00/frame; 18038.128 us/call
```

Wall min/avg/max: 61.755/62.388/63.101 ms. Window mean render is 44.351 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 17–32)

```text
frame                                     62.421 ms 37.453 Mcyc 100.0%
  pov_preserve_half                        0.147 ms  0.088 Mcyc   0.2%
    x1.00/frame; 147.420 us/call
  gx_draw_particles                        1.027 ms  0.616 Mcyc   1.6%
    x1.00/frame; 1026.952 us/call
  gx_particle_step                         0.972 ms  0.583 Mcyc   1.6%
    x1.00/frame; 972.493 us/call
  gx_timeline_step                         0.087 ms  0.052 Mcyc   0.1%
    x1.00/frame; 87.311 us/call
  canvas_clear                             0.086 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.725 us/call
  canvas_buffer_wait                      60.074 ms 36.044 Mcyc  96.2%
    x1.00/frame; 60073.891 us/call
```

Wall min/avg/max: 62.197/62.421/62.489 ms. Window mean render is 2.347 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.06/frame
  min/avg/max 0.538/1.776/30.343 us; cpu 3.27%
  isr_pack         144.00/frame
    min/avg/max 6.043/7.633/10.658 us; cpu 1.75%
  isr_dma_submit   144.00/frame
    min/avg/max 0.650/0.940/1.031 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999994 us
isr_wake       n=18433 cyc min/avg/max=323/1066/18206 ns min/avg/max=538/1776/30343 total=32754 us cpu=3.27%
isr_pack       n=2304 cyc min/avg/max=3626/4580/6395 ns min/avg/max=6043/7633/10658 total=17591 us cpu=1.75%
isr_dma_submit n=2304 cyc min/avg/max=390/564/619 ns min/avg/max=650/940/1031 total=2168 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.27% leaves approximately 60.456 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `gx_particle_step` — 27.101 ms/frame, 43.52% of this root window.
2. `gx_draw_particles` — 16.454 ms/frame, 26.42% of this root window.
3. `pov_preserve_half` — 0.145 ms/frame, 0.23% of this root window.
4. `gx_timeline_step` — 0.087 ms/frame, 0.14% of this root window.
5. `canvas_clear` — 0.086 ms/frame, 0.14% of this root window.

README cells: peak 🟢 45.126, spilled 🟢 0/927 (0.00%).

README setup frame 1: 1.482 ms. Captured 2026-11-04 22:44 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Galaxies  windows=58
  [PASS] captured >=3 windows (58)
  [PASS] every window names one effect and resolution (Galaxies 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 58 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 58 windows are short)
  [PASS] scope 'gx_particle_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.5 ppm, frames 865-880)

VALID: Galaxies
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Galaxies`, `HS_PROFILE_WINDOW=16`. `just profile Galaxies` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
