# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/gsreactiondiffusion_ship.log`. Captured 2026-11-04 22:59 America/Los_Angeles on COM4. Replaces profile_gsreactiondiffusion_teensy_2026-10-30.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GSReactionDiffusion, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh GSReactionDiffusion profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:84784, data:249048, headers:9208   free for files:1688576
teensy_size:    RAM1: variables:314272, code:17592, padding:15176   free for local variables:177248
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `921deabc0eedcea5155e6febea778ca7045a4884f0df1b0ffc8117d5df62e149`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 385–416, root cyc/600 versus measured wall sum differs by **0.098 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   GSReactionDiffusion      worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 10,318 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **48.568 ms** (frame 395), **0/447 spilled** (0.00%). Setup frame 1: 24.927 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 40.991 ms/frame. `grd_render` averages 40.663 ms/frame over complete windows; its worst window is 47.611 ms/frame (frames 385–416).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 20.201–48.568 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     60.176 ms 36.106 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.082 Mcyc   0.2%
    x0.97/frame; 140.504 us/call
  grd_render                              24.705 ms 14.823 Mcyc  41.1%
    grd_shader_draw                        9.778 ms  5.867 Mcyc  16.2%
      x6.00/frame; 1629.674 us/call
    grd_cull_flags                         0.717 ms  0.430 Mcyc   1.2%
      x6.00/frame; 119.564 us/call
    grd_publish                            1.962 ms  1.177 Mcyc   3.3%
      x6.00/frame; 327.028 us/call
    grd_physics                            5.160 ms  3.096 Mcyc   8.6%
      x6.00/frame; 859.994 us/call
    grd_pigment                            3.074 ms  1.844 Mcyc   5.1%
      x1.00/frame; 3073.751 us/call
    rd_lattice_basis                       0.113 ms  0.068 Mcyc   0.2%
      x1.00/frame; 113.380 us/call
    grd_color_noise                        2.992 ms  1.795 Mcyc   5.0%
      x1.00/frame; 2991.726 us/call
  rd_timeline_step                         0.045 ms  0.027 Mcyc   0.1%
    x1.00/frame; 45.487 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.2%
    x1.00/frame; 91.989 us/call
  canvas_buffer_wait                      35.177 ms 21.106 Mcyc  58.5%
    x1.00/frame; 35177.088 us/call
grd_color_palette                          0.430 ms  0.258 Mcyc   0.7%
  x0.03/frame; 13746.010 us/call
grd_seed_reaction                          0.487 ms  0.292 Mcyc   0.8%
  x0.03/frame; 15569.295 us/call
```

Wall min/avg/max: 21.446/60.176/64.898 ms. Window mean render is 24.999 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 385–416)

```text
frame                                     62.422 ms 37.453 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.081 Mcyc   0.2%
    x1.00/frame; 135.708 us/call
  grd_render                              47.611 ms 28.567 Mcyc  76.3%
    grd_shader_draw                       28.016 ms 16.810 Mcyc  44.9%
      x6.00/frame; 4669.402 us/call
    grd_cull_flags                         0.740 ms  0.444 Mcyc   1.2%
      x6.00/frame; 123.275 us/call
    grd_publish                            4.672 ms  2.803 Mcyc   7.5%
      x6.00/frame; 778.659 us/call
    grd_physics                            5.184 ms  3.110 Mcyc   8.3%
      x6.00/frame; 863.970 us/call
    grd_pigment                            4.232 ms  2.539 Mcyc   6.8%
      x1.00/frame; 4232.397 us/call
    rd_lattice_basis                       0.111 ms  0.066 Mcyc   0.2%
      x1.00/frame; 110.539 us/call
    grd_color_noise                        3.015 ms  1.809 Mcyc   4.8%
      x1.00/frame; 3014.718 us/call
  rd_timeline_step                         0.045 ms  0.027 Mcyc   0.1%
    x1.00/frame; 44.698 us/call
  canvas_clear                             0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.148 us/call
  canvas_buffer_wait                      14.520 ms  8.712 Mcyc  23.3%
    x1.00/frame; 14520.482 us/call
```

Wall min/avg/max: 61.801/62.422/63.184 ms. Window mean render is 47.902 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.568 ms 37.541 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 136.823 us/call
  grd_render                              30.280 ms 18.168 Mcyc  48.4%
    grd_shader_draw                       14.853 ms  8.912 Mcyc  23.7%
      x6.00/frame; 2475.552 us/call
    grd_cull_flags                         0.715 ms  0.429 Mcyc   1.1%
      x6.00/frame; 119.175 us/call
    grd_publish                            1.960 ms  1.176 Mcyc   3.1%
      x6.00/frame; 326.620 us/call
    grd_physics                            5.159 ms  3.095 Mcyc   8.2%
      x6.00/frame; 859.857 us/call
    grd_pigment                            3.578 ms  2.147 Mcyc   5.7%
      x1.00/frame; 3578.234 us/call
    rd_lattice_basis                       0.110 ms  0.066 Mcyc   0.2%
      x1.00/frame; 110.262 us/call
    grd_color_noise                        2.999 ms  1.800 Mcyc   4.8%
      x1.00/frame; 2999.376 us/call
  rd_timeline_step                         0.044 ms  0.026 Mcyc   0.1%
    x1.00/frame; 43.821 us/call
  canvas_clear                             0.089 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.927 us/call
  canvas_buffer_wait                      31.999 ms 19.199 Mcyc  51.1%
    x1.00/frame; 31998.536 us/call
```

Wall min/avg/max: 56.848/62.568/68.403 ms. Window mean render is 30.570 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1147.84/frame
  min/avg/max 0.538/1.880/24.483 us; cpu 3.46%
  isr_pack         143.47/frame
    min/avg/max 6.338/8.525/16.388 us; cpu 1.96%
  isr_dma_submit   143.47/frame
    min/avg/max 0.691/0.946/4.565 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1992659 us
isr_wake       n=36731 cyc min/avg/max=323/1128/14690 ns min/avg/max=538/1880/24483 total=69076 us cpu=3.46%
isr_pack       n=4591 cyc min/avg/max=3803/5115/9833 ns min/avg/max=6338/8525/16388 total=39141 us cpu=1.96%
isr_dma_submit n=4591 cyc min/avg/max=415/568/2739 ns min/avg/max=691/946/4565 total=4348 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.46% leaves approximately 60.338 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `grd_render` — 47.611 ms/frame, 76.27% of this root window.
2. `grd_shader_draw` — 28.016 ms/frame, 44.88% of this root window.
3. `grd_physics` — 5.184 ms/frame, 8.30% of this root window.
4. `grd_publish` — 4.672 ms/frame, 7.48% of this root window.
5. `grd_pigment` — 4.232 ms/frame, 6.78% of this root window.

README cells: peak 🟢 48.568, spilled 🟢 0/447 (0.00%).

README setup frame 1: 24.927 ms. Captured 2026-11-04 22:59 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=GSReactionDiffusion  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (GSReactionDiffusion 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'grd_render' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.1 ppm, frames 385-416)

VALID: GSReactionDiffusion
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`. `just profile GSReactionDiffusion` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
