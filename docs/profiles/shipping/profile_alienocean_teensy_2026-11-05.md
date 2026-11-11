# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienOcean`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/alienocean_ship.log`. Captured 2026-11-04 23:49 America/Los_Angeles on COM4. Replaces profile_alienocean_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienOcean, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh AlienOcean profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:100832, data:152616, headers:8696   free for files:1769472
teensy_size:    RAM1: variables:314144, code:28904, padding:3864   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `0f130e4542c5781aed28adaa71ca582e5237411062aeaff7922eefb30633ac55`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **1.561 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   AlienOcean               worst 2 lines/set (0 sets over)  gated 20,474 B  overflow  1,812 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **26.435 ms** (frame 433), **0/447 spilled** (0.00%). Setup frame 1: 48.819 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 24.774 ms/frame. `fx_shader_draw` averages 22.416 ms/frame over complete windows; its worst window is 23.058 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 23.220–26.435 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     60.750 ms 36.450 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x0.97/frame; 140.916 us/call
  fx_shader_draw                          23.058 ms 13.835 Mcyc  38.0%
    fx_prepare_column                      0.795 ms  0.477 Mcyc   1.3%
      x144.00/frame; 5.519 us/call
  fx_palette_base_bake                     0.189 ms  0.114 Mcyc   0.3%
    fx_palette_rgb_bake                    0.059 ms  0.036 Mcyc   0.1%
      x1.00/frame; 59.296 us/call
  fx_prepare_frame                         0.401 ms  0.240 Mcyc   0.7%
    x1.00/frame; 400.573 us/call
  fx_advance                               1.189 ms  0.713 Mcyc   2.0%
    palette_cycle_step                     1.102 ms  0.661 Mcyc   1.8%
      x1.00/frame; 1101.689 us/call
  fx_timeline_step                         0.093 ms  0.056 Mcyc   0.2%
    x1.00/frame; 92.889 us/call
  canvas_clear                             0.101 ms  0.061 Mcyc   0.2%
    x1.00/frame; 101.058 us/call
  canvas_buffer_wait                      35.343 ms 21.206 Mcyc  58.2%
    x1.00/frame; 35342.649 us/call
```

Wall min/avg/max: 23.456/60.750/65.304 ms. Window mean render is 25.408 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 417–448)

```text
frame                                     62.396 ms 37.437 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.091 us/call
  fx_shader_draw                          22.292 ms 13.375 Mcyc  35.7%
    fx_prepare_column                      0.796 ms  0.478 Mcyc   1.3%
      x144.00/frame; 5.528 us/call
  fx_palette_base_bake                     0.515 ms  0.309 Mcyc   0.8%
    fx_palette_rgb_bake                    0.386 ms  0.232 Mcyc   0.6%
      x1.00/frame; 386.103 us/call
  fx_prepare_frame                         0.328 ms  0.197 Mcyc   0.5%
    x1.00/frame; 328.001 us/call
  fx_advance                               1.181 ms  0.708 Mcyc   1.9%
    palette_cycle_step                     1.155 ms  0.693 Mcyc   1.9%
      x1.00/frame; 1154.896 us/call
  fx_timeline_step                         0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.924 us/call
  canvas_clear                             0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 96.821 us/call
  canvas_buffer_wait                      37.523 ms 22.514 Mcyc  60.1%
    x1.00/frame; 37523.391 us/call
```

Wall min/avg/max: 59.591/62.395/65.195 ms. Window mean render is 24.872 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 65–96)

```text
frame                                     62.376 ms 37.425 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.304 us/call
  fx_shader_draw                          22.364 ms 13.418 Mcyc  35.9%
    fx_prepare_column                      0.793 ms  0.476 Mcyc   1.3%
      x144.00/frame; 5.506 us/call
  fx_palette_base_bake                     0.258 ms  0.155 Mcyc   0.4%
    fx_palette_rgb_bake                    0.128 ms  0.077 Mcyc   0.2%
      x1.00/frame; 128.431 us/call
  fx_prepare_frame                         0.312 ms  0.187 Mcyc   0.5%
    x1.00/frame; 312.278 us/call
  fx_advance                               1.172 ms  0.703 Mcyc   1.9%
    palette_cycle_step                     1.145 ms  0.687 Mcyc   1.8%
      x1.00/frame; 1145.067 us/call
  fx_timeline_step                         0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.034 us/call
  canvas_clear                             0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 96.907 us/call
  canvas_buffer_wait                      37.712 ms 22.627 Mcyc  60.5%
    x1.00/frame; 37711.578 us/call
```

Wall min/avg/max: 59.235/62.375/65.577 ms. Window mean render is 24.664 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.22/frame
  min/avg/max 0.495/1.966/21.835 us; cpu 3.62%
  isr_pack         144.03/frame
    min/avg/max 7.835/9.456/11.995 us; cpu 2.17%
  isr_dma_submit   144.03/frame
    min/avg/max 0.778/0.938/1.508 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000254 us
isr_wake       n=36871 cyc min/avg/max=297/1180/13101 ns min/avg/max=495/1966/21835 total=72548 us cpu=3.62%
isr_pack       n=4609 cyc min/avg/max=4701/5674/7197 ns min/avg/max=7835/9456/11995 total=43589 us cpu=2.17%
isr_dma_submit n=4609 cyc min/avg/max=467/563/905 ns min/avg/max=778/938/1508 total=4328 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.62% leaves approximately 60.237 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 22.292 ms/frame, 35.73% of this root window.
2. `fx_advance` — 1.181 ms/frame, 1.89% of this root window.
3. `palette_cycle_step` — 1.155 ms/frame, 1.85% of this root window.
4. `fx_prepare_column` — 0.796 ms/frame, 1.28% of this root window.
5. `fx_palette_base_bake` — 0.515 ms/frame, 0.83% of this root window.

README cells: peak 🟢 26.435, spilled 🟢 0/447 (0.00%).

README setup frame 1: 48.819 ms. Captured 2026-11-04 23:49 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=AlienOcean  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (AlienOcean 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.6 ppm, frames 1-32)

VALID: AlienOcean
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienOcean`, `HS_PROFILE_WINDOW=32`. `just profile AlienOcean` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
