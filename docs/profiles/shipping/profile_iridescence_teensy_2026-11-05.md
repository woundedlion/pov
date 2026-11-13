# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Iridescence`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/iridescence_ship.log`. Captured 2026-11-05 00:22 America/Los_Angeles on COM4. Replaces profile_iridescence_teensy_2026-11-04.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Iridescence, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Iridescence profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:110680, data:157228, headers:8572   free for files:1755136
teensy_size:    RAM1: variables:314144, code:28600, padding:4168   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `b86df7d761002785c6955d7cbd23843b54ed4a22e1a1fd7ac47dbb7e3484e9e7`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **1.273 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Iridescence              worst 2 lines/set (0 sets over)  gated 20,476 B  overflow  7,360 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **49.094 ms** (frame 25), **0/927 spilled** (0.00%). Setup frame 1: 91.161 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 47.498 ms/frame. `fx_shader_draw` averages 42.948 ms/frame over complete windows; its worst window is 44.128 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 45.814–49.094 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.781 ms 37.669 Mcyc 100.0%
  pov_preserve_half                        0.134 ms  0.080 Mcyc   0.2%
    x0.97/frame; 137.898 us/call
  fx_shader_draw                          44.128 ms 26.477 Mcyc  70.3%
    fx_prepare_column                      1.594 ms  0.957 Mcyc   2.5%
      x144.00/frame; 11.072 us/call
  fx_palette_base_bake                     0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.308 us/call
  fx_prepare_frame                         3.158 ms  1.895 Mcyc   5.0%
    x1.00/frame; 3158.243 us/call
  fx_advance                               1.102 ms  0.661 Mcyc   1.8%
    palette_cycle_step                     1.015 ms  0.609 Mcyc   1.6%
      x1.00/frame; 1015.342 us/call
  fx_timeline_step                         0.018 ms  0.011 Mcyc   0.0%
    x1.00/frame; 17.820 us/call
  canvas_clear                             0.105 ms  0.063 Mcyc   0.2%
    x1.00/frame; 105.144 us/call
  canvas_buffer_wait                      13.762 ms  8.257 Mcyc  21.9%
    x1.00/frame; 13762.477 us/call
```

Wall min/avg/max: 46.439/62.781/91.163 ms. Window mean render is 49.019 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 449–480)

```text
frame                                     62.371 ms 37.423 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 139.011 us/call
  fx_shader_draw                          42.806 ms 25.684 Mcyc  68.6%
    fx_prepare_column                      1.623 ms  0.974 Mcyc   2.6%
      x144.00/frame; 11.268 us/call
  fx_palette_base_bake                     0.089 ms  0.053 Mcyc   0.1%
    x1.00/frame; 88.552 us/call
  fx_prepare_frame                         3.018 ms  1.811 Mcyc   4.8%
    x1.00/frame; 3017.740 us/call
  fx_advance                               0.922 ms  0.553 Mcyc   1.5%
    palette_cycle_step                     0.900 ms  0.540 Mcyc   1.4%
      x1.00/frame; 900.257 us/call
  fx_timeline_step                         0.018 ms  0.011 Mcyc   0.0%
    x1.00/frame; 17.827 us/call
  canvas_clear                             0.101 ms  0.061 Mcyc   0.2%
    x1.00/frame; 100.976 us/call
  canvas_buffer_wait                      14.994 ms  8.996 Mcyc  24.0%
    x1.00/frame; 14993.777 us/call
```

Wall min/avg/max: 60.309/62.371/64.463 ms. Window mean render is 47.377 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.97/frame
  min/avg/max 0.506/1.990/18.615 us; cpu 3.66%
  isr_pack         144.00/frame
    min/avg/max 8.623/9.323/11.831 us; cpu 2.14%
  isr_dma_submit   144.00/frame
    min/avg/max 0.606/0.936/6.926 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999845 us
isr_wake       n=36863 cyc min/avg/max=304/1194/11169 ns min/avg/max=506/1990/18615 total=73387 us cpu=3.66%
isr_pack       n=4608 cyc min/avg/max=5174/5594/7099 ns min/avg/max=8623/9323/11831 total=42967 us cpu=2.14%
isr_dma_submit n=4608 cyc min/avg/max=364/562/4156 ns min/avg/max=606/936/6926 total=4318 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.66% leaves approximately 60.213 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 44.128 ms/frame, 70.29% of this root window.
2. `fx_prepare_frame` — 3.158 ms/frame, 5.03% of this root window.
3. `fx_prepare_column` — 1.594 ms/frame, 2.54% of this root window.
4. `fx_advance` — 1.102 ms/frame, 1.75% of this root window.
5. `palette_cycle_step` — 1.015 ms/frame, 1.62% of this root window.

README cells: peak 🟢 49.094, spilled 🟢 0/927 (0.00%).

README setup frame 1: 91.161 ms. Captured 2026-11-05 00:22 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Iridescence  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (Iridescence 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.3 ppm, frames 1-32)

VALID: Iridescence
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Iridescence`, `HS_PROFILE_WINDOW=32`. `just profile Iridescence` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
