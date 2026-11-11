# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienCore`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/aliencore_ship.log`. Captured 2026-11-04 23:50 America/Los_Angeles on COM4. Replaces profile_aliencore_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienCore, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh AlienCore profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:101464, data:152616, headers:9088   free for files:1768448
teensy_size:    RAM1: variables:314144, code:28456, padding:4312   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `a92c952457bf9ae915cb6e6f8b7501b0f1530074975993cec649020a7c9dee61`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **0.496 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   AlienCore                worst 2 lines/set (0 sets over)  gated 20,478 B  overflow  3,048 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **25.172 ms** (frame 429), **0/447 spilled** (0.00%). Setup frame 1: 40.847 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 20.656 ms/frame. `fx_shader_draw` averages 18.113 ms/frame over complete windows; its worst window is 18.572 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 18.785–25.172 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     60.371 ms 36.222 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x0.97/frame; 141.468 us/call
  fx_shader_draw                          18.572 ms 11.143 Mcyc  30.8%
    fx_prepare_column                      0.651 ms  0.391 Mcyc   1.1%
      x144.00/frame; 4.520 us/call
  fx_palette_base_bake                     0.193 ms  0.116 Mcyc   0.3%
    fx_palette_rgb_bake                    0.077 ms  0.046 Mcyc   0.1%
      x1.00/frame; 76.791 us/call
  fx_prepare_frame                         0.421 ms  0.252 Mcyc   0.7%
    x1.00/frame; 420.696 us/call
  fx_advance                               1.141 ms  0.684 Mcyc   1.9%
    palette_cycle_step                     1.045 ms  0.627 Mcyc   1.7%
      x1.00/frame; 1044.550 us/call
  fx_timeline_step                         0.076 ms  0.046 Mcyc   0.1%
    x1.00/frame; 76.012 us/call
  canvas_clear                             0.096 ms  0.058 Mcyc   0.2%
    x1.00/frame; 96.006 us/call
  canvas_buffer_wait                      39.499 ms 23.699 Mcyc  65.4%
    x1.00/frame; 39499.164 us/call
```

Wall min/avg/max: 18.920/60.370/65.616 ms. Window mean render is 20.872 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 417–448)

```text
frame                                     62.397 ms 37.438 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 144.673 us/call
  fx_shader_draw                          18.145 ms 10.887 Mcyc  29.1%
    fx_prepare_column                      0.648 ms  0.389 Mcyc   1.0%
      x144.00/frame; 4.502 us/call
  fx_palette_base_bake                     0.797 ms  0.478 Mcyc   1.3%
    fx_palette_rgb_bake                    0.681 ms  0.409 Mcyc   1.1%
      x1.00/frame; 680.924 us/call
  fx_prepare_frame                         0.382 ms  0.229 Mcyc   0.6%
    x1.00/frame; 381.549 us/call
  fx_advance                               1.153 ms  0.692 Mcyc   1.8%
    palette_cycle_step                     1.118 ms  0.671 Mcyc   1.8%
      x1.00/frame; 1117.711 us/call
  fx_timeline_step                         0.078 ms  0.047 Mcyc   0.1%
    x1.00/frame; 77.715 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.536 us/call
  canvas_buffer_wait                      41.368 ms 24.821 Mcyc  66.3%
    x1.00/frame; 41367.829 us/call
```

Wall min/avg/max: 56.696/62.397/68.120 ms. Window mean render is 21.029 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.399 ms 37.439 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 142.033 us/call
  fx_shader_draw                          18.063 ms 10.838 Mcyc  28.9%
    fx_prepare_column                      0.652 ms  0.391 Mcyc   1.0%
      x144.00/frame; 4.525 us/call
  fx_palette_base_bake                     0.255 ms  0.153 Mcyc   0.4%
    fx_palette_rgb_bake                    0.138 ms  0.083 Mcyc   0.2%
      x1.00/frame; 137.865 us/call
  fx_prepare_frame                         0.331 ms  0.199 Mcyc   0.5%
    x1.00/frame; 331.281 us/call
  fx_advance                               1.136 ms  0.681 Mcyc   1.8%
    palette_cycle_step                     1.100 ms  0.660 Mcyc   1.8%
      x1.00/frame; 1100.179 us/call
  fx_timeline_step                         0.076 ms  0.045 Mcyc   0.1%
    x1.00/frame; 75.563 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 92.149 us/call
  canvas_buffer_wait                      42.074 ms 25.244 Mcyc  67.4%
    x1.00/frame; 42073.604 us/call
```

Wall min/avg/max: 59.098/62.399/65.687 ms. Window mean render is 20.325 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.09/frame
  min/avg/max 0.533/1.986/28.588 us; cpu 3.66%
  isr_pack         144.00/frame
    min/avg/max 8.653/9.418/11.946 us; cpu 2.17%
  isr_dma_submit   144.00/frame
    min/avg/max 0.678/0.943/7.300 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2000066 us
isr_wake       n=36867 cyc min/avg/max=320/1192/17153 ns min/avg/max=533/1986/28588 total=73271 us cpu=3.66%
isr_pack       n=4608 cyc min/avg/max=5192/5651/7168 ns min/avg/max=8653/9418/11946 total=43406 us cpu=2.17%
isr_dma_submit n=4608 cyc min/avg/max=407/566/4380 ns min/avg/max=678/943/7300 total=4352 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.66% leaves approximately 60.213 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 18.145 ms/frame, 29.08% of this root window.
2. `fx_advance` — 1.153 ms/frame, 1.85% of this root window.
3. `palette_cycle_step` — 1.118 ms/frame, 1.79% of this root window.
4. `fx_palette_base_bake` — 0.797 ms/frame, 1.28% of this root window.
5. `fx_palette_rgb_bake` — 0.681 ms/frame, 1.09% of this root window.

README cells: peak 🟢 25.172, spilled 🟢 0/447 (0.00%).

README setup frame 1: 40.847 ms. Captured 2026-11-04 23:50 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=AlienCore  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (AlienCore 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.5 ppm, frames 1-32)

VALID: AlienCore
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienCore`, `HS_PROFILE_WINDOW=32`. `just profile AlienCore` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
