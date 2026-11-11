# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Fishbowl`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/fishbowl_ship.log`. Captured 2026-11-04 22:52 America/Los_Angeles on COM4. Replaces profile_fishbowl_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Fishbowl, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Fishbowl profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:140280, data:149688, headers:9040   free for files:1732608
teensy_size:    RAM1: variables:314048, code:38568, padding:26968   free for local variables:144704
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `8aa0b709567f2e80c321272d50725d369da48ba6999bf236d509c294882cdd05`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 737–768, root cyc/600 versus measured wall sum differs by **2.591 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Fishbowl                 worst 2 lines/set (0 sets over)  gated 20,474 B  overflow 48,472 B  full-set worst 5
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **28.932 ms** (frame 215), **0/927 spilled** (0.00%). Setup frame 1: 1.755 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 22.044 ms/frame. `fish_build_vertices` averages 11.883 ms/frame over complete windows; its worst window is 12.681 ms/frame (frames 737–768).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 1.495–28.932 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     58.777 ms 35.266 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x0.97/frame; 148.067 us/call
  fish_multiline_draw                      2.075 ms  1.245 Mcyc   3.5%
    x1.00/frame; 2075.309 us/call
  fish_build_vertices                      1.865 ms  1.119 Mcyc   3.2%
    x1.00/frame; 1865.383 us/call
  fish_noise_prepare                       0.002 ms  0.001 Mcyc   0.0%
    x1.00/frame; 2.435 us/call
  fish_timeline_step                       0.170 ms  0.102 Mcyc   0.3%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.178 us/call
  canvas_clear                             0.096 ms  0.057 Mcyc   0.2%
    x1.00/frame; 95.687 us/call
  canvas_buffer_wait                      54.271 ms 32.563 Mcyc  92.3%
    x1.00/frame; 54271.450 us/call
```

Wall min/avg/max: 1.514/58.777/63.901 ms. Window mean render is 4.506 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 193–224)

```text
frame                                     62.375 ms 37.425 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.149 us/call
  fish_multiline_draw                     11.274 ms  6.764 Mcyc  18.1%
    x1.00/frame; 11273.780 us/call
  fish_build_vertices                     12.627 ms  7.576 Mcyc  20.2%
    x1.00/frame; 12627.196 us/call
  fish_noise_prepare                       0.002 ms  0.001 Mcyc   0.0%
    x1.00/frame; 2.099 us/call
  fish_timeline_step                       0.187 ms  0.112 Mcyc   0.3%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.177 us/call
  canvas_clear                             0.091 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.392 us/call
  canvas_buffer_wait                      37.385 ms 22.431 Mcyc  59.9%
    x1.00/frame; 37385.372 us/call
```

Wall min/avg/max: 58.075/62.374/67.030 ms. Window mean render is 24.990 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 737–768)

```text
frame                                     62.373 ms 37.424 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.442 us/call
  fish_multiline_draw                     10.162 ms  6.097 Mcyc  16.3%
    x1.00/frame; 10162.278 us/call
  fish_build_vertices                     12.681 ms  7.609 Mcyc  20.3%
    x1.00/frame; 12680.950 us/call
  fish_noise_prepare                       0.003 ms  0.002 Mcyc   0.0%
    x1.00/frame; 2.790 us/call
  fish_timeline_step                       0.187 ms  0.112 Mcyc   0.3%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.177 us/call
  canvas_clear                             0.091 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.749 us/call
  canvas_buffer_wait                      38.445 ms 23.067 Mcyc  61.6%
    x1.00/frame; 38445.148 us/call
```

Wall min/avg/max: 59.676/62.372/64.330 ms. Window mean render is 23.928 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.490 ms 37.494 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.400 us/call
  fish_multiline_draw                      4.274 ms  2.565 Mcyc   6.8%
    x1.00/frame; 4274.243 us/call
  fish_build_vertices                      5.376 ms  3.226 Mcyc   8.6%
    x1.00/frame; 5376.340 us/call
  fish_noise_prepare                       0.002 ms  0.001 Mcyc   0.0%
    x1.00/frame; 2.151 us/call
  fish_timeline_step                       0.181 ms  0.109 Mcyc   0.3%
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.177 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.453 us/call
  canvas_buffer_wait                      52.099 ms 31.260 Mcyc  83.4%
    x1.00/frame; 52099.424 us/call
```

Wall min/avg/max: 61.015/62.490/63.910 ms. Window mean render is 10.391 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.50/frame
  min/avg/max 0.558/1.683/25.350 us; cpu 3.10%
  isr_pack         143.94/frame
    min/avg/max 5.991/6.700/9.568 us; cpu 1.54%
  isr_dma_submit   143.94/frame
    min/avg/max 0.681/0.940/9.346 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1998986 us
isr_wake       n=36848 cyc min/avg/max=335/1010/15210 ns min/avg/max=558/1683/25350 total=62078 us cpu=3.10%
isr_pack       n=4606 cyc min/avg/max=3595/4020/5741 ns min/avg/max=5991/6700/9568 total=30867 us cpu=1.54%
isr_dma_submit n=4606 cyc min/avg/max=409/564/5608 ns min/avg/max=681/940/9346 total=4333 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.10% leaves approximately 60.562 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fish_build_vertices` — 12.627 ms/frame, 20.24% of this root window.
2. `fish_multiline_draw` — 11.274 ms/frame, 18.07% of this root window.
3. `fish_timeline_step` — 0.187 ms/frame, 0.30% of this root window.
4. `pov_preserve_half` — 0.141 ms/frame, 0.23% of this root window.
5. `canvas_clear` — 0.091 ms/frame, 0.15% of this root window.

README cells: peak 🟢 28.932, spilled 🟢 0/927 (0.00%).

README setup frame 1: 1.755 ms. Captured 2026-11-04 22:52 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Fishbowl  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (Fishbowl 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fish_build_vertices' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (2.6 ppm, frames 737-768)

VALID: Fishbowl
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Fishbowl`, `HS_PROFILE_WINDOW=32`. `just profile Fishbowl` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
