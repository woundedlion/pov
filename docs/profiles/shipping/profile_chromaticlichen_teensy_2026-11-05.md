# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ChromaticLichen`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/chromaticlichen_ship.log`. Captured 2026-11-05 00:05 America/Los_Angeles on COM4. Replaces profile_chromaticlichen_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ChromaticLichen, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh ChromaticLichen profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:102488, data:151736, headers:8944   free for files:1768448
teensy_size:    RAM1: variables:314144, code:29960, padding:2808   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `872a0c64b650db314de3f929fc7c9bef7be4ee6d14bd511e1113a553103919bd`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **0.480 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   ChromaticLichen          worst 2 lines/set (0 sets over)  gated 20,478 B  overflow  2,936 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **40.138 ms** (frame 257), **0/447 spilled** (0.00%). Setup frame 1: 63.665 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 35.543 ms/frame. `fx_shader_draw` averages 29.603 ms/frame over complete windows; its worst window is 30.204 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 32.763–40.138 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     61.540 ms 36.924 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x0.97/frame; 141.621 us/call
  fx_shader_draw                          30.204 ms 18.122 Mcyc  49.1%
    fx_prepare_column                      0.638 ms  0.383 Mcyc   1.0%
      x144.00/frame; 4.430 us/call
  fx_palette_base_bake                     0.229 ms  0.137 Mcyc   0.4%
    fx_palette_rgb_bake                    0.109 ms  0.065 Mcyc   0.2%
      x1.00/frame; 108.650 us/call
  fx_prepare_frame                         3.285 ms  1.971 Mcyc   5.3%
    x1.00/frame; 3285.188 us/call
  fx_advance                               1.427 ms  0.856 Mcyc   2.3%
    palette_cycle_step                     1.313 ms  0.788 Mcyc   2.1%
      x1.00/frame; 1313.438 us/call
  fx_timeline_step                         0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 87.618 us/call
  canvas_clear                             0.098 ms  0.059 Mcyc   0.2%
    x1.00/frame; 97.830 us/call
  canvas_buffer_wait                      25.867 ms 15.520 Mcyc  42.0%
    x1.00/frame; 25867.499 us/call
```

Wall min/avg/max: 32.987/61.539/66.855 ms. Window mean render is 35.672 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 257–288)

```text
frame                                     62.410 ms 37.446 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.607 us/call
  fx_shader_draw                          29.447 ms 17.668 Mcyc  47.2%
    fx_prepare_column                      0.636 ms  0.382 Mcyc   1.0%
      x144.00/frame; 4.418 us/call
  fx_palette_base_bake                     1.065 ms  0.639 Mcyc   1.7%
    fx_palette_rgb_bake                    0.942 ms  0.565 Mcyc   1.5%
      x1.00/frame; 942.362 us/call
  fx_prepare_frame                         3.356 ms  2.014 Mcyc   5.4%
    x1.00/frame; 3355.894 us/call
  fx_advance                               1.212 ms  0.727 Mcyc   1.9%
    palette_cycle_step                     1.182 ms  0.709 Mcyc   1.9%
      x1.00/frame; 1182.105 us/call
  fx_timeline_step                         0.084 ms  0.051 Mcyc   0.1%
    x1.00/frame; 84.440 us/call
  canvas_clear                             0.093 ms  0.056 Mcyc   0.1%
    x1.00/frame; 93.342 us/call
  canvas_buffer_wait                      26.819 ms 16.091 Mcyc  43.0%
    x1.00/frame; 26818.897 us/call
```

Wall min/avg/max: 55.984/62.409/68.874 ms. Window mean render is 35.591 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 417–448)

```text
frame                                     62.407 ms 37.444 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.709 us/call
  fx_shader_draw                          29.273 ms 17.564 Mcyc  46.9%
    fx_prepare_column                      0.635 ms  0.381 Mcyc   1.0%
      x144.00/frame; 4.411 us/call
  fx_palette_base_bake                     0.684 ms  0.410 Mcyc   1.1%
    fx_palette_rgb_bake                    0.561 ms  0.337 Mcyc   0.9%
      x1.00/frame; 561.108 us/call
  fx_prepare_frame                         3.244 ms  1.946 Mcyc   5.2%
    x1.00/frame; 3244.117 us/call
  fx_advance                               1.145 ms  0.687 Mcyc   1.8%
    palette_cycle_step                     1.117 ms  0.670 Mcyc   1.8%
      x1.00/frame; 1117.272 us/call
  fx_timeline_step                         0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.015 us/call
  canvas_clear                             0.094 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.412 us/call
  canvas_buffer_wait                      27.548 ms 16.529 Mcyc  44.1%
    x1.00/frame; 27547.699 us/call
```

Wall min/avg/max: 59.387/62.407/65.361 ms. Window mean render is 34.860 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.06/frame
  min/avg/max 0.498/1.958/22.310 us; cpu 3.61%
  isr_pack         144.00/frame
    min/avg/max 7.950/9.126/11.620 us; cpu 2.10%
  isr_dma_submit   144.00/frame
    min/avg/max 0.663/0.931/1.915 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999939 us
isr_wake       n=36866 cyc min/avg/max=299/1175/13386 ns min/avg/max=498/1958/22310 total=72221 us cpu=3.61%
isr_pack       n=4608 cyc min/avg/max=4770/5476/6972 ns min/avg/max=7950/9126/11620 total=42058 us cpu=2.10%
isr_dma_submit n=4608 cyc min/avg/max=398/559/1149 ns min/avg/max=663/931/1915 total=4299 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.61% leaves approximately 60.244 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 29.447 ms/frame, 47.18% of this root window.
2. `fx_prepare_frame` — 3.356 ms/frame, 5.38% of this root window.
3. `fx_advance` — 1.212 ms/frame, 1.94% of this root window.
4. `palette_cycle_step` — 1.182 ms/frame, 1.89% of this root window.
5. `fx_palette_base_bake` — 1.065 ms/frame, 1.71% of this root window.

README cells: peak 🟢 40.138, spilled 🟢 0/447 (0.00%).

README setup frame 1: 63.665 ms. Captured 2026-11-05 00:05 America/Los_Angeles.

Pullback shader: also affected by the 1a3244da8..dc2954193 cached-phase/color fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=ChromaticLichen  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (ChromaticLichen 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.5 ppm, frames 1-32)

VALID: ChromaticLichen
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`. `just profile ChromaticLichen` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
