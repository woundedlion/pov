# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile CosmicEyeball`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/cosmiceyeball_ship.log`. Captured 2026-11-05 00:40 America/Los_Angeles on COM4. Replaces profile_cosmiceyeball_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | CosmicEyeball, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh CosmicEyeball profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:93488, data:152152, headers:8312   free for files:1777664
teensy_size:    RAM1: variables:314144, code:25896, padding:6872   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `4702554d383864d990847972ad7ff703f52684bb39a44943fe6673553af34b55`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **1.247 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   CosmicEyeball            worst 2 lines/set (0 sets over)  gated 18,098 B  overflow      0 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **25.370 ms** (frame 539), **0/927 spilled** (0.00%). Setup frame 1: 43.803 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 23.488 ms/frame. `fx_shader_draw` averages 21.297 ms/frame over complete windows; its worst window is 21.884 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 21.793–25.370 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     60.552 ms 36.331 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x0.97/frame; 142.291 us/call
  fx_shader_draw                          21.884 ms 13.130 Mcyc  36.1%
    fx_prepare_column                      0.825 ms  0.495 Mcyc   1.4%
      x144.00/frame; 5.731 us/call
  fx_palette_base_bake                     0.149 ms  0.089 Mcyc   0.2%
    fx_palette_rgb_bake                    0.046 ms  0.028 Mcyc   0.1%
      x1.00/frame; 46.269 us/call
  fx_prepare_frame                         0.297 ms  0.178 Mcyc   0.5%
    x1.00/frame; 296.709 us/call
  fx_advance                               1.118 ms  0.671 Mcyc   1.8%
    palette_cycle_step                     1.036 ms  0.622 Mcyc   1.7%
      x1.00/frame; 1036.251 us/call
  fx_timeline_step                         0.077 ms  0.046 Mcyc   0.1%
    x1.00/frame; 76.760 us/call
  canvas_clear                             0.091 ms  0.055 Mcyc   0.2%
    x1.00/frame; 91.108 us/call
  canvas_buffer_wait                      36.586 ms 21.951 Mcyc  60.4%
    x1.00/frame; 36585.575 us/call
```

Wall min/avg/max: 21.978/60.552/65.582 ms. Window mean render is 23.967 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 513–544)

```text
frame                                     62.368 ms 37.421 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 142.903 us/call
  fx_shader_draw                          21.305 ms 12.783 Mcyc  34.2%
    fx_prepare_column                      0.800 ms  0.480 Mcyc   1.3%
      x144.00/frame; 5.553 us/call
  fx_palette_base_bake                     0.259 ms  0.156 Mcyc   0.4%
    fx_palette_rgb_bake                    0.148 ms  0.089 Mcyc   0.2%
      x1.00/frame; 148.442 us/call
  fx_prepare_frame                         0.296 ms  0.178 Mcyc   0.5%
    x1.00/frame; 296.229 us/call
  fx_advance                               1.090 ms  0.654 Mcyc   1.7%
    palette_cycle_step                     1.072 ms  0.643 Mcyc   1.7%
      x1.00/frame; 1071.874 us/call
  fx_timeline_step                         0.074 ms  0.044 Mcyc   0.1%
    x1.00/frame; 73.962 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.363 us/call
  canvas_buffer_wait                      38.916 ms 23.350 Mcyc  62.4%
    x1.00/frame; 38915.926 us/call
```

Wall min/avg/max: 58.992/62.368/65.768 ms. Window mean render is 23.452 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 609–640)

```text
frame                                     62.398 ms 37.439 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 142.776 us/call
  fx_shader_draw                          21.253 ms 12.752 Mcyc  34.1%
    fx_prepare_column                      0.821 ms  0.493 Mcyc   1.3%
      x144.00/frame; 5.702 us/call
  fx_palette_base_bake                     0.169 ms  0.102 Mcyc   0.3%
    fx_palette_rgb_bake                    0.061 ms  0.037 Mcyc   0.1%
      x1.00/frame; 61.493 us/call
  fx_prepare_frame                         0.297 ms  0.178 Mcyc   0.5%
    x1.00/frame; 297.199 us/call
  fx_advance                               1.065 ms  0.639 Mcyc   1.7%
    palette_cycle_step                     1.047 ms  0.628 Mcyc   1.7%
      x1.00/frame; 1046.609 us/call
  fx_timeline_step                         0.074 ms  0.044 Mcyc   0.1%
    x1.00/frame; 74.084 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.215 us/call
  canvas_buffer_wait                      39.106 ms 23.464 Mcyc  62.7%
    x1.00/frame; 39106.109 us/call
```

Wall min/avg/max: 59.051/62.398/65.728 ms. Window mean render is 23.292 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.06/frame
  min/avg/max 0.496/1.983/26.968 us; cpu 3.65%
  isr_pack         144.00/frame
    min/avg/max 8.771/9.551/12.381 us; cpu 2.20%
  isr_dma_submit   144.00/frame
    min/avg/max 0.733/0.933/4.861 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999977 us
isr_wake       n=36866 cyc min/avg/max=298/1190/16181 ns min/avg/max=496/1983/26968 total=73159 us cpu=3.65%
isr_pack       n=4608 cyc min/avg/max=5263/5731/7429 ns min/avg/max=8771/9551/12381 total=44021 us cpu=2.20%
isr_dma_submit n=4608 cyc min/avg/max=440/560/2917 ns min/avg/max=733/933/4861 total=4303 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.65% leaves approximately 60.219 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 21.305 ms/frame, 34.16% of this root window.
2. `fx_advance` — 1.090 ms/frame, 1.75% of this root window.
3. `palette_cycle_step` — 1.072 ms/frame, 1.72% of this root window.
4. `fx_prepare_column` — 0.800 ms/frame, 1.28% of this root window.
5. `fx_prepare_frame` — 0.296 ms/frame, 0.47% of this root window.

README cells: peak 🟢 25.370, spilled 🟢 0/927 (0.00%).

README setup frame 1: 43.803 ms. Captured 2026-11-05 00:40 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=CosmicEyeball  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (CosmicEyeball 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.2 ppm, frames 1-32)

VALID: CosmicEyeball
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=CosmicEyeball`, `HS_PROFILE_WINDOW=32`. `just profile CosmicEyeball` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
