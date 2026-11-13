# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Millefiori`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/millefiori_ship.log`. Captured 2026-11-04 23:58 America/Los_Angeles on COM4. Replaces profile_millefiori_teensy_2026-10-30.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Millefiori, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Millefiori profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:106704, data:153944, headers:8664   free for files:1762304
teensy_size:    RAM1: variables:314112, code:33192, padding:32344   free for local variables:144640
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `bdef2d46aabb3233d417229d5e70b41a4da97779f9cd6cb0f0c5d19c03ae341a`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 193–224, root cyc/600 versus measured wall sum differs by **0.092 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Millefiori               worst 2 lines/set (0 sets over)  gated 20,476 B  overflow  2,912 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **51.758 ms** (frame 389), **0/927 spilled** (0.00%). Setup frame 1: 78.463 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 46.247 ms/frame. `fx_shader_draw` averages 43.852 ms/frame over complete windows; its worst window is 46.890 ms/frame (frames 193–224).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 38.280–51.758 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.317 ms 37.390 Mcyc 100.0%
  pov_preserve_half                        0.134 ms  0.080 Mcyc   0.2%
    x0.97/frame; 138.209 us/call
  fx_shader_draw                          42.644 ms 25.587 Mcyc  68.4%
    fx_prepare_column                      1.412 ms  0.847 Mcyc   2.3%
      x144.00/frame; 9.802 us/call
  fx_palette_base_bake                     0.115 ms  0.069 Mcyc   0.2%
    x1.00/frame; 114.835 us/call
  fx_prepare_frame                         0.579 ms  0.347 Mcyc   0.9%
    x1.00/frame; 578.717 us/call
  fx_advance                               1.166 ms  0.700 Mcyc   1.9%
    palette_cycle_step                     1.062 ms  0.637 Mcyc   1.7%
      x1.00/frame; 1061.630 us/call
  fx_timeline_step                         0.104 ms  0.063 Mcyc   0.2%
    x1.00/frame; 104.366 us/call
  canvas_clear                             0.101 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.823 us/call
  canvas_buffer_wait                      17.236 ms 10.342 Mcyc  27.7%
    x1.00/frame; 17236.372 us/call
```

Wall min/avg/max: 39.855/62.317/78.467 ms. Window mean render is 45.081 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 385–416)

```text
frame                                     62.455 ms 37.473 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.653 us/call
  fx_shader_draw                          44.004 ms 26.403 Mcyc  70.5%
    fx_prepare_column                      1.406 ms  0.844 Mcyc   2.3%
      x144.00/frame; 9.764 us/call
  fx_palette_base_bake                     0.116 ms  0.069 Mcyc   0.2%
    x1.00/frame; 115.651 us/call
  fx_prepare_frame                         0.548 ms  0.329 Mcyc   0.9%
    x1.00/frame; 547.865 us/call
  fx_advance                               1.165 ms  0.699 Mcyc   1.9%
    palette_cycle_step                     1.128 ms  0.677 Mcyc   1.8%
      x1.00/frame; 1128.289 us/call
  fx_timeline_step                         0.108 ms  0.065 Mcyc   0.2%
    x1.00/frame; 107.631 us/call
  canvas_clear                             0.098 ms  0.059 Mcyc   0.2%
    x1.00/frame; 97.948 us/call
  canvas_buffer_wait                      16.041 ms  9.625 Mcyc  25.7%
    x1.00/frame; 16041.139 us/call
```

Wall min/avg/max: 53.867/62.454/70.482 ms. Window mean render is 46.414 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 193–224)

```text
frame                                     62.406 ms 37.444 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 139.115 us/call
  fx_shader_draw                          46.890 ms 28.134 Mcyc  75.1%
    fx_prepare_column                      1.422 ms  0.853 Mcyc   2.3%
      x144.00/frame; 9.873 us/call
  fx_palette_base_bake                     0.112 ms  0.067 Mcyc   0.2%
    x1.00/frame; 112.498 us/call
  fx_prepare_frame                         0.719 ms  0.431 Mcyc   1.2%
    x1.00/frame; 718.938 us/call
  fx_advance                               1.110 ms  0.666 Mcyc   1.8%
    palette_cycle_step                     1.070 ms  0.642 Mcyc   1.7%
      x1.00/frame; 1070.070 us/call
  fx_timeline_step                         0.109 ms  0.065 Mcyc   0.2%
    x1.00/frame; 109.101 us/call
  canvas_clear                             0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 95.401 us/call
  canvas_buffer_wait                      12.992 ms  7.795 Mcyc  20.8%
    x1.00/frame; 12991.817 us/call
```

Wall min/avg/max: 60.009/62.406/64.823 ms. Window mean render is 49.415 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 737–768)

```text
frame                                     62.324 ms 37.394 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.792 us/call
  fx_shader_draw                          38.459 ms 23.075 Mcyc  61.7%
    fx_prepare_column                      1.410 ms  0.846 Mcyc   2.3%
      x144.00/frame; 9.793 us/call
  fx_palette_base_bake                     0.116 ms  0.069 Mcyc   0.2%
    x1.00/frame; 115.712 us/call
  fx_prepare_frame                         0.495 ms  0.297 Mcyc   0.8%
    x1.00/frame; 494.919 us/call
  fx_advance                               1.147 ms  0.688 Mcyc   1.8%
    palette_cycle_step                     1.109 ms  0.665 Mcyc   1.8%
      x1.00/frame; 1108.663 us/call
  fx_timeline_step                         0.105 ms  0.063 Mcyc   0.2%
    x1.00/frame; 105.315 us/call
  canvas_clear                             0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 95.127 us/call
  canvas_buffer_wait                      21.527 ms 12.916 Mcyc  34.5%
    x1.00/frame; 21526.664 us/call
```

Wall min/avg/max: 58.426/62.323/66.154 ms. Window mean render is 40.797 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1154.59/frame
  min/avg/max 0.585/1.991/30.043 us; cpu 3.67%
  isr_pack         144.31/frame
    min/avg/max 8.631/9.336/12.221 us; cpu 2.15%
  isr_dma_submit   144.31/frame
    min/avg/max 0.665/0.936/1.045 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2004342 us
isr_wake       n=36947 cyc min/avg/max=351/1195/18026 ns min/avg/max=585/1991/30043 total=73612 us cpu=3.67%
isr_pack       n=4618 cyc min/avg/max=5179/5602/7333 ns min/avg/max=8631/9336/12221 total=43120 us cpu=2.15%
isr_dma_submit n=4618 cyc min/avg/max=399/562/627 ns min/avg/max=665/936/1045 total=4333 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.67% leaves approximately 60.206 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 44.004 ms/frame, 70.46% of this root window.
2. `fx_prepare_column` — 1.406 ms/frame, 2.25% of this root window.
3. `fx_advance` — 1.165 ms/frame, 1.86% of this root window.
4. `palette_cycle_step` — 1.128 ms/frame, 1.81% of this root window.
5. `fx_prepare_frame` — 0.548 ms/frame, 0.88% of this root window.

README cells: peak 🟢 51.758, spilled 🟢 0/927 (0.00%).

README setup frame 1: 78.463 ms. Captured 2026-11-04 23:58 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Millefiori  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (Millefiori 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.1 ppm, frames 193-224)

VALID: Millefiori
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Millefiori`, `HS_PROFILE_WINDOW=32`. `just profile Millefiori` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
