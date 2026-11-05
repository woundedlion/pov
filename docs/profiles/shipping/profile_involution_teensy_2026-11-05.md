# Involution on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Involution`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/involution_ship.log`. Captured 2026-11-05 00:15 America/Los_Angeles on COM4. Replaces profile_involution_teensy_2026-11-04.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Involution, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Involution profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:112032, data:159544, headers:9000   free for files:1751040
teensy_size:    RAM1: variables:314144, code:29448, padding:3320   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `d8be8b78eddf14f92415527f8abbf03a398787ed80972b670d107c0e47ead66d`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **2.545 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Involution               worst 2 lines/set (0 sets over)  gated 20,472 B  overflow  1,732 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **47.347 ms** (frame 848), **0/927 spilled** (0.00%). Setup frame 1: 83.449 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 45.482 ms/frame. `fx_shader_draw` averages 43.789 ms/frame over complete windows; its worst window is 44.606 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 43.468–47.347 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.487 ms 37.492 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.97/frame; 139.266 us/call
  fx_shader_draw                          44.606 ms 26.764 Mcyc  71.4%
    fx_prepare_column                      3.186 ms  1.911 Mcyc   5.1%
      x144.00/frame; 22.123 us/call
  fx_palette_base_bake                     0.083 ms  0.050 Mcyc   0.1%
    x1.00/frame; 83.049 us/call
  fx_prepare_frame                         0.113 ms  0.068 Mcyc   0.2%
    x1.00/frame; 113.265 us/call
  fx_advance                               1.073 ms  0.644 Mcyc   1.7%
    palette_cycle_step                     0.972 ms  0.583 Mcyc   1.6%
      x1.00/frame; 972.490 us/call
  fx_timeline_step                         0.068 ms  0.041 Mcyc   0.1%
    x1.00/frame; 67.703 us/call
  canvas_clear                             0.094 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.171 us/call
  canvas_buffer_wait                      16.030 ms  9.618 Mcyc  25.7%
    x1.00/frame; 16030.081 us/call
```

Wall min/avg/max: 43.606/62.486/83.453 ms. Window mean render is 46.457 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 833–864)

```text
frame                                     62.392 ms 37.435 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 139.069 us/call
  fx_shader_draw                          43.991 ms 26.395 Mcyc  70.5%
    fx_prepare_column                      3.730 ms  2.238 Mcyc   6.0%
      x144.00/frame; 25.904 us/call
  fx_palette_base_bake                     0.082 ms  0.049 Mcyc   0.1%
    x1.00/frame; 82.251 us/call
  fx_prepare_frame                         0.027 ms  0.016 Mcyc   0.0%
    x1.00/frame; 27.021 us/call
  fx_advance                               1.063 ms  0.638 Mcyc   1.7%
    palette_cycle_step                     1.025 ms  0.615 Mcyc   1.6%
      x1.00/frame; 1024.543 us/call
  fx_timeline_step                         0.067 ms  0.040 Mcyc   0.1%
    x1.00/frame; 66.707 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.937 us/call
  canvas_buffer_wait                      16.652 ms  9.991 Mcyc  26.7%
    x1.00/frame; 16652.245 us/call
```

Wall min/avg/max: 59.698/62.391/65.207 ms. Window mean render is 45.740 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.405 ms 37.443 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.370 us/call
  fx_shader_draw                          43.301 ms 25.981 Mcyc  69.4%
    fx_prepare_column                      3.165 ms  1.899 Mcyc   5.1%
      x144.00/frame; 21.976 us/call
  fx_palette_base_bake                     0.083 ms  0.050 Mcyc   0.1%
    x1.00/frame; 82.961 us/call
  fx_prepare_frame                         0.029 ms  0.017 Mcyc   0.0%
    x1.00/frame; 28.849 us/call
  fx_advance                               1.064 ms  0.638 Mcyc   1.7%
    palette_cycle_step                     1.025 ms  0.615 Mcyc   1.6%
      x1.00/frame; 1025.266 us/call
  fx_timeline_step                         0.067 ms  0.040 Mcyc   0.1%
    x1.00/frame; 67.269 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.027 us/call
  canvas_buffer_wait                      17.353 ms 10.412 Mcyc  27.8%
    x1.00/frame; 17353.326 us/call
```

Wall min/avg/max: 59.705/62.404/64.990 ms. Window mean render is 45.051 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.00/frame
  min/avg/max 0.630/2.003/28.688 us; cpu 3.69%
  isr_pack         144.00/frame
    min/avg/max 7.453/9.415/12.066 us; cpu 2.16%
  isr_dma_submit   144.00/frame
    min/avg/max 0.633/0.935/1.030 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999859 us
isr_wake       n=36864 cyc min/avg/max=378/1202/17213 ns min/avg/max=630/2003/28688 total=73891 us cpu=3.69%
isr_pack       n=4608 cyc min/avg/max=4472/5649/7240 ns min/avg/max=7453/9415/12066 total=43390 us cpu=2.16%
isr_dma_submit n=4608 cyc min/avg/max=380/561/618 ns min/avg/max=633/935/1030 total=4312 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.69% leaves approximately 60.194 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 43.991 ms/frame, 70.51% of this root window.
2. `fx_prepare_column` — 3.730 ms/frame, 5.98% of this root window.
3. `fx_advance` — 1.063 ms/frame, 1.70% of this root window.
4. `palette_cycle_step` — 1.025 ms/frame, 1.64% of this root window.
5. `pov_preserve_half` — 0.139 ms/frame, 0.22% of this root window.

README cells: peak 🟢 47.347, spilled 🟢 0/927 (0.00%).

README setup frame 1: 83.449 ms. Captured 2026-11-05 00:15 America/Los_Angeles.

Pullback shader: also affected by the 1a3244da8..dc2954193 cached-phase/color fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Involution  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (Involution 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (2.5 ppm, frames 1-32)

VALID: Involution
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Involution`, `HS_PROFILE_WINDOW=32`. `just profile Involution` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
