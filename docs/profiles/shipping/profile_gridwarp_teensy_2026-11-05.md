# GridWarp on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GridWarp`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/gridwarp_ship.log`. Captured 2026-11-05 00:19 America/Los_Angeles on COM4. Replaces profile_gridwarp_teensy_2026-11-04.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridWarp, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh GridWarp profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:107856, data:154864, headers:8640   free for files:1760256
teensy_size:    RAM1: variables:314144, code:29672, padding:3096   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `dc524b85f28cc1f1a177c95af646db87f2fc82789d7d7e0dabe5b1dcbe75aa53`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 353–384, root cyc/600 versus measured wall sum differs by **0.209 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   GridWarp                 worst 2 lines/set (0 sets over)  gated 20,480 B  overflow  7,168 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **49.281 ms** (frame 245), **0/927 spilled** (0.00%). Setup frame 1: 74.267 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 44.050 ms/frame. `fx_shader_draw` averages 38.505 ms/frame over complete windows; its worst window is 39.195 ms/frame (frames 353–384).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 41.618–49.281 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.102 ms 37.261 Mcyc 100.0%
  pov_preserve_half                        0.134 ms  0.080 Mcyc   0.2%
    x0.97/frame; 138.013 us/call
  fx_shader_draw                          39.193 ms 23.516 Mcyc  63.1%
    fx_prepare_column                      4.432 ms  2.659 Mcyc   7.1%
      x144.00/frame; 30.781 us/call
  fx_palette_base_bake                     0.260 ms  0.156 Mcyc   0.4%
    fx_palette_rgb_bake                    0.120 ms  0.072 Mcyc   0.2%
      x1.00/frame; 120.234 us/call
  fx_prepare_frame                         3.380 ms  2.028 Mcyc   5.4%
    x1.00/frame; 3379.568 us/call
  fx_advance                               0.983 ms  0.590 Mcyc   1.6%
    palette_cycle_step                     0.895 ms  0.537 Mcyc   1.4%
      x1.00/frame; 895.099 us/call
  fx_timeline_step                         0.104 ms  0.062 Mcyc   0.2%
    x1.00/frame; 104.030 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.513 us/call
  canvas_buffer_wait                      17.730 ms 10.638 Mcyc  28.5%
    x1.00/frame; 17730.013 us/call
```

Wall min/avg/max: 41.954/62.101/74.269 ms. Window mean render is 44.372 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 225–256)

```text
frame                                     62.362 ms 37.417 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 137.624 us/call
  fx_shader_draw                          38.658 ms 23.195 Mcyc  62.0%
    fx_prepare_column                      4.460 ms  2.676 Mcyc   7.2%
      x144.00/frame; 30.971 us/call
  fx_palette_base_bake                     1.337 ms  0.802 Mcyc   2.1%
    fx_palette_rgb_bake                    1.192 ms  0.715 Mcyc   1.9%
      x1.00/frame; 1192.359 us/call
  fx_prepare_frame                         3.442 ms  2.065 Mcyc   5.5%
    x1.00/frame; 3442.198 us/call
  fx_advance                               0.921 ms  0.552 Mcyc   1.5%
    palette_cycle_step                     0.886 ms  0.531 Mcyc   1.4%
      x1.00/frame; 885.591 us/call
  fx_timeline_step                         0.109 ms  0.065 Mcyc   0.2%
    x1.00/frame; 109.013 us/call
  canvas_clear                             0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 87.547 us/call
  canvas_buffer_wait                      17.449 ms 10.469 Mcyc  28.0%
    x1.00/frame; 17449.035 us/call
```

Wall min/avg/max: 55.908/62.362/68.836 ms. Window mean render is 44.913 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 353–384)

```text
frame                                     62.362 ms 37.417 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.693 us/call
  fx_shader_draw                          39.195 ms 23.517 Mcyc  62.9%
    fx_prepare_column                      4.462 ms  2.677 Mcyc   7.2%
      x144.00/frame; 30.984 us/call
  fx_palette_base_bake                     1.496 ms  0.898 Mcyc   2.4%
    fx_palette_rgb_bake                    1.351 ms  0.810 Mcyc   2.2%
      x1.00/frame; 1350.774 us/call
  fx_prepare_frame                         3.377 ms  2.026 Mcyc   5.4%
    x1.00/frame; 3376.520 us/call
  fx_advance                               0.986 ms  0.592 Mcyc   1.6%
    palette_cycle_step                     0.952 ms  0.571 Mcyc   1.5%
      x1.00/frame; 952.109 us/call
  fx_timeline_step                         0.105 ms  0.063 Mcyc   0.2%
    x1.00/frame; 105.411 us/call
  canvas_clear                             0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.314 us/call
  canvas_buffer_wait                      16.760 ms 10.056 Mcyc  26.9%
    x1.00/frame; 16759.898 us/call
```

Wall min/avg/max: 56.350/62.362/68.452 ms. Window mean render is 45.602 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 577–608)

```text
frame                                     62.436 ms 37.461 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 137.316 us/call
  fx_shader_draw                          38.360 ms 23.016 Mcyc  61.4%
    fx_prepare_column                      4.461 ms  2.676 Mcyc   7.1%
      x144.00/frame; 30.977 us/call
  fx_palette_base_bake                     0.198 ms  0.119 Mcyc   0.3%
    fx_palette_rgb_bake                    0.063 ms  0.038 Mcyc   0.1%
      x1.00/frame; 62.862 us/call
  fx_prepare_frame                         3.292 ms  1.975 Mcyc   5.3%
    x1.00/frame; 3291.891 us/call
  fx_advance                               0.923 ms  0.554 Mcyc   1.5%
    palette_cycle_step                     0.889 ms  0.534 Mcyc   1.4%
      x1.00/frame; 889.195 us/call
  fx_timeline_step                         0.106 ms  0.064 Mcyc   0.2%
    x1.00/frame; 106.077 us/call
  canvas_clear                             0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 87.614 us/call
  canvas_buffer_wait                      19.106 ms 11.463 Mcyc  30.6%
    x1.00/frame; 19105.746 us/call
```

Wall min/avg/max: 59.688/62.435/65.054 ms. Window mean render is 43.330 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.38/frame
  min/avg/max 0.538/2.015/30.401 us; cpu 3.71%
  isr_pack         144.03/frame
    min/avg/max 8.770/9.490/12.213 us; cpu 2.18%
  isr_dma_submit   144.03/frame
    min/avg/max 0.643/0.940/5.310 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000500 us
isr_wake       n=36876 cyc min/avg/max=323/1209/18241 ns min/avg/max=538/2015/30401 total=74349 us cpu=3.71%
isr_pack       n=4609 cyc min/avg/max=5262/5694/7328 ns min/avg/max=8770/9490/12213 total=43741 us cpu=2.18%
isr_dma_submit n=4609 cyc min/avg/max=386/564/3186 ns min/avg/max=643/940/5310 total=4337 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.71% leaves approximately 60.181 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 38.658 ms/frame, 61.99% of this root window.
2. `fx_prepare_column` — 4.460 ms/frame, 7.15% of this root window.
3. `fx_prepare_frame` — 3.442 ms/frame, 5.52% of this root window.
4. `fx_palette_base_bake` — 1.337 ms/frame, 2.14% of this root window.
5. `fx_palette_rgb_bake` — 1.192 ms/frame, 1.91% of this root window.

README cells: peak 🟢 49.281, spilled 🟢 0/927 (0.00%).

README setup frame 1: 74.267 ms. Captured 2026-11-05 00:19 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=GridWarp  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (GridWarp 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.2 ppm, frames 353-384)

VALID: GridWarp
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridWarp`, `HS_PROFILE_WINDOW=32`. `just profile GridWarp` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
