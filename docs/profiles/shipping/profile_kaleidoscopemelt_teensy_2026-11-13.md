# KaleidoscopeMelt on-device profile — Teensy 4.0, segmented mode (2026-11-13, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeMelt`).
Raw capture: `build/prof/kaleidoscopemelt_ship41b.log` in `C:/work/temp/Holosphere-20261112-ashcloud-itcm`. Captured 2026-11-13 07:06 America/Los_Angeles on COM4 (serial 6362390). First shipping-roster report. Replaces profile_kaleidoscopemelt_selective_o3_teensy_2026-10-31.md (excluded set).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4 serial 6362390, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: shipping `-Os` base, newlib-nano, `-ffast-math -fno-finite-math-only`; selective `HS_O3` column draw in `core/render/scan/shader.h` and pullback shader/pipeline paths |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeMelt 288×144, single-entry playlist, source `6ab24361ce9699afe9b9365c0c7b6f360abde12d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 64 s capture; default one-hour epoch, no dwell compression |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261112-ashcloud-itcm HS_DEVICE_WAIT=1800 HS_TEENSY_PORT=COM4 HS_PROFILE_OUT=build/prof/kaleidoscopemelt_ship41b.log bash tools/profile_one.sh KaleidoscopeMelt profile 64 32` |

Image size (measured single-entry Profile image):

```text
teensy_size:   FLASH: code:102456, data:152004, headers:8708   free for files:1768448
teensy_size:    RAM1: variables:314112, code:35096, padding:30440   free for local variables:144640
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Shipping Phantasm image (41 entries):

```text
teensy_size:   FLASH: code:1183144, data:792876, headers:8492   free for files:47104
teensy_size:    RAM1: variables:313760, code:195656, padding:952   free for local variables:13920
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

The archived Phantasm build passes teensy-gate and hot-layout: RAM1 code is 195,656 B of 196,608 B, with 952 B remaining in six ITCM banks. RAM1 variables are 313,760 B; free local-variable space is 13,920 B against the 12,288 B stack floor.

Profile ELF SHA-256: `f5cea20762d3d8704a87d9e5ae055941a7448fedcf9707f83bd14b63d967b51a`. Phantasm ELF SHA-256: `050da474f94d77520d8f320508fee88272d55c44764393466064e3158204a162`. Both archived ELFs in `build/prof/artifacts/kaleidoscopemelt_ship41b_6ab24361ce96_f5cea20762d3` were rehashed for this report. Adjacent `kaleidoscopemelt_ship41b.provenance`, `kaleidoscopemelt_ship41b_envdump.txt` and `kaleidoscopemelt_ship41b_phantasm_envdump.txt` retain source, flags and compiler/ABI: Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 20251203.

Exactness cross-check: window frames 193–224, root 1,197,458,164 cyc ÷ 600 MHz = 1995763.607 us versus measured wall sum 1,995,768 us, within **2.201 ppm** (`python tools/parse_profile.py build/prof/kaleidoscopemelt_ship41b.log validate`; 5 ppm limit).

Hot-layout line in the measured Profile build (PASS):

```text
  ok   KaleidoscopeMelt         worst 2 lines/set (0 sets over)  gated 19,294 B  overflow      0 B  full-set worst 2
```

Corresponding line in the 41-entry Phantasm build (PASS):

```text
  ok   KaleidoscopeMelt         worst 2 lines/set (0 sets over)  gated 20,178 B  overflow      0 B  full-set worst 2
```

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer, setup included): peak frame render 77.194 ms, spilled 1/992 frames (0.10%). `fx_shader_draw` averages 41.832 ms/f; its worst window is 43.730 ms/f (frames 193–224).

**Steady render, frame 1 excluded:** peak **50.140 ms** (frame 932), **0/991 spilled** (0.00%). Setup frame 1: 77.194 ms. Complete windows cover frames 1–992; steady statistics use frames 2–992. Trailing rows without a complete counter window are excluded; subsequent transitions remain included. Steady mean render is 46.452 ms/frame.

A display window is 62.5 ms at 480 RPM; the effect renders one quadrant ≈ 10,368 px. All measured steady frames fit 16 fps, with 12.360 ms peak margin; steady render spans 42.423–50.140 ms. Setup exceeds one display window and is reported separately. `canvas_buffer_wait` is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the default shipping timeline runs through startup and a continuously evolving shader; no preset markers occur. The windows below show startup, the steady peak, the worst shader window and the lowest-cost steady window where distinct. Counter values are per-frame cycle-derived averages; percentages use the root frame, and nested costs overlap.

### Startup (window frames 1–32)

```text
frame                          62.269 ms  37.361 Mcyc 100.0%
  pov_preserve_half           133.997 us  80.398 kcyc   0.2%
    x0.97/frame; 138.319 us/call
  fx_shader_draw               41.337 ms  24.802 Mcyc  66.4%
    fx_prepare_column         762.911 us 457.746 kcyc   1.2%
      x144.00/frame; 5.298 us/call
  fx_palette_base_bake        327.145 us 196.287 kcyc   0.5%
    x1.00/frame; 327.145 us/call
  fx_prepare_frame              2.926 ms   1.756 Mcyc   4.7%
    x1.00/frame; 2926.273 us/call
  fx_advance                    1.047 ms 628.105 kcyc   1.7%
    palette_cycle_step          1.016 ms 609.399 kcyc   1.6%
      x1.00/frame; 1015.666 us/call
  fx_timeline_step             95.060 us  57.036 kcyc   0.2%
    x1.00/frame; 95.060 us/call
  canvas_clear                 96.808 us  58.085 kcyc   0.2%
    x1.00/frame; 96.808 us/call
  canvas_buffer_wait           16.188 ms   9.713 Mcyc  26.0%
    x1.00/frame; 16187.713 us/call
```

Wall min/avg/max = 43.942/62.269/77.197 ms. Window mean render is 46.081 ms. This window includes setup frame 1; its maximum is excluded from the steady ranking. Later windows retain the changing shader workload without that setup draw.

### Peak-render regime (window frames 929–960)

```text
frame                          62.357 ms  37.414 Mcyc 100.0%
  pov_preserve_half           136.658 us  81.995 kcyc   0.2%
    x1.00/frame; 136.658 us/call
  fx_shader_draw               42.651 ms  25.591 Mcyc  68.4%
    fx_prepare_column         774.869 us 464.922 kcyc   1.2%
      x144.00/frame; 5.381 us/call
  fx_palette_base_bake        338.881 us 203.329 kcyc   0.5%
    x1.00/frame; 338.881 us/call
  fx_prepare_frame              2.982 ms   1.789 Mcyc   4.8%
    x1.00/frame; 2982.202 us/call
  fx_advance                  904.919 us 542.951 kcyc   1.5%
    palette_cycle_step        874.977 us 524.986 kcyc   1.4%
      x1.00/frame; 874.977 us/call
  fx_timeline_step             93.558 us  56.135 kcyc   0.2%
    x1.00/frame; 93.558 us/call
  canvas_clear                 92.157 us  55.294 kcyc   0.1%
    x1.00/frame; 92.157 us/call
  canvas_buffer_wait           15.046 ms   9.028 Mcyc  24.1%
    x1.00/frame; 15046.218 us/call
```

Wall min/avg/max = 58.189/62.357/66.744 ms. Window mean render is 47.311 ms. This window contains the worst steady frame; the counter means do not equal that single-frame peak. The shader dominates live work, while the buffer wait fills the remaining display interval.

### Worst dominant-scope regime (window frames 193–224)

```text
frame                          62.368 ms  37.421 Mcyc 100.0%
  pov_preserve_half           138.300 us  82.980 kcyc   0.2%
    x1.00/frame; 138.300 us/call
  fx_shader_draw               43.730 ms  26.238 Mcyc  70.1%
    fx_prepare_column         775.080 us 465.048 kcyc   1.2%
      x144.00/frame; 5.382 us/call
  fx_palette_base_bake        336.162 us 201.697 kcyc   0.5%
    x1.00/frame; 336.162 us/call
  fx_prepare_frame              2.937 ms   1.762 Mcyc   4.7%
    x1.00/frame; 2936.780 us/call
  fx_advance                  889.214 us 533.529 kcyc   1.4%
    palette_cycle_step        858.648 us 515.189 kcyc   1.4%
      x1.00/frame; 858.648 us/call
  fx_timeline_step             94.255 us  56.553 kcyc   0.2%
    x1.00/frame; 94.255 us/call
  canvas_clear                 92.693 us  55.616 kcyc   0.1%
    x1.00/frame; 92.693 us/call
  canvas_buffer_wait           14.046 ms   8.427 Mcyc  22.5%
    x1.00/frame; 14045.501 us/call
```

Wall min/avg/max = 60.323/62.368/64.235 ms. Window mean render is 48.322 ms. This window has the largest mean shader cost, although its single-frame render peak is lower than the peak-render window. The shader and preparation scopes explain the reduced idle interval relative to the lowest-cost regime.

### Lowest-cost steady regime (window frames 705–736)

```text
frame                          62.376 ms  37.425 Mcyc 100.0%
  pov_preserve_half           140.237 us  84.142 kcyc   0.2%
    x1.00/frame; 140.237 us/call
  fx_shader_draw               39.747 ms  23.848 Mcyc  63.7%
    fx_prepare_column         764.372 us 458.623 kcyc   1.2%
      x144.00/frame; 5.308 us/call
  fx_palette_base_bake        335.760 us 201.456 kcyc   0.5%
    x1.00/frame; 335.760 us/call
  fx_prepare_frame              2.967 ms   1.780 Mcyc   4.8%
    x1.00/frame; 2967.356 us/call
  fx_advance                  913.531 us 548.119 kcyc   1.5%
    palette_cycle_step        882.365 us 529.419 kcyc   1.4%
      x1.00/frame; 882.365 us/call
  fx_timeline_step             90.469 us  54.281 kcyc   0.1%
    x1.00/frame; 90.469 us/call
  canvas_clear                 92.385 us  55.431 kcyc   0.1%
    x1.00/frame; 92.385 us/call
  canvas_buffer_wait           17.972 ms  10.783 Mcyc  28.8%
    x1.00/frame; 17972.041 us/call
```

Wall min/avg/max = 59.058/62.376/65.715 ms. Window mean render is 44.404 ms. This is the lowest mean-render window after startup. Its larger buffer wait reflects the lower live workload at the same 16 fps cadence.

### Per-pixel figures

Per-pixel blend counters are disabled in this shipping timing image. Blended px/frame, quadrant coverage, cycles/blend and scan cycles/blended pixel are unavailable; geometric quadrant size is not substituted for measured blend counts.

## Column-ISR / DMA marshaling cost

Last complete window (frames 961–992); separately measured ISR interval = 2,000,405 us. Rates divide counts by 32 frames; CPU shares use that ISR interval, which includes report output.

```text
isr_wake        1152.31/frame
  min/avg/max 0.285/1.758/22.228 us  cpu 3.24%
  isr_pack        144.03/frame
    min/avg/max 8.090/8.735/11.348 us  cpu 2.01%
  isr_dma_submit  144.03/frame
    min/avg/max 0.601/0.915/2.466 us  cpu 0.21%
```

- CPU submission averages 0.915 us/call versus 8.735 us/call for packing (0.21% versus 2.01% CPU).
- Submission measures CPU work; asynchronous wire-transfer time is not measured.
- Inclusive `isr_wake` CPU share is 3.24%, leaving approximately 60.475 ms of foreground CPU per 62.5 ms window. Pack and submit run inside wake; use wake alone. DMA-completion and other interrupts are not measured. Render scopes already absorb interrupts, so no second ISR subtraction applies to their 62.5 ms budget. Every measured steady regime requires 1.000× speedup (already fits); setup alone would require 1.235× to fit one window.

## Summary ranking

1. `fx_shader_draw` — 68.40% of the frame, 42.651 ms: column shader draw.
2. `fx_prepare_frame` — 4.78% of the frame, 2.982 ms: frame preparation.
3. `fx_advance` — 1.45% of the frame, 0.905 ms: animation and palette advancement.
4. `palette_cycle_step` — 1.40% of the frame, 0.875 ms: palette cycling nested in advancement.
5. `fx_prepare_column` — 1.24% of the frame, 0.775 ms: column preparation nested in shader draw.

Ranked costs use the peak-render window, frames 929–960; nested scopes are inclusive. Native/WASM results are not a matched hardware timing baseline; this capture supplies no paired native/WASM measurement.

README cells: peak 🟢 50.140, spilled 🟢 0/991 (0.00%).

README setup frame 1: 77.194 ms. Captured 2026-11-13 07:06 America/Los_Angeles.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` first-caller parenting can hide its subtree when that parent has no calls; calls approximate blended pixels when enabled. Per-pixel scope overhead applies when enabled; per-pixel counters are disabled here.
- Shipping selective `HS_O3` column-draw and pullback shader/pipeline paths are active; this is the selective-O3 shipping configuration.
- No dwell-compression, ordered-cycle or transition-speed overrides were used; the default one-hour epoch exceeds this 64 s capture. These knobs control coverage rather than per-frame cost.
- Provenance records source `6ab24361ce9699afe9b9365c0c7b6f360abde12d`; working-tree cleanliness is not recorded in the run artifacts.
- This segment-0 master capture bounds the measured default-control interval; other segment roles and exhaustive parameter combinations are not measured.

Parser validation:

```text
effect=KaleidoscopeMelt  windows=31
  [PASS] captured >=3 windows (31)
  [PASS] every window names one effect and resolution (KaleidoscopeMelt 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 31 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 31 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (2.2 ppm, frames 193-224)

VALID: KaleidoscopeMelt
```

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeMelt`, `HS_PROFILE_WINDOW=32`; `just profile KaleidoscopeMelt` = build + flash + capture through the shared-device lock. The Setup command reproduces this capture's duration, port and output path.
