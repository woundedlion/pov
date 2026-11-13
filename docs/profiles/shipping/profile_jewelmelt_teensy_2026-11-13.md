# JewelMelt on-device profile — Teensy 4.0, segmented mode (2026-11-13, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile JewelMelt`).
Raw capture: `build/prof/jewelmelt_ship41b.log` in `C:/work/temp/Holosphere-20261112-ashcloud-itcm`. Captured 2026-11-13 07:02 America/Los_Angeles on COM4 (serial 6362390). First shipping-roster report. Replaces `excluded/profile_jewelmelt_selective_o3_teensy_2026-11-03.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4 serial 6362390, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: shipping `-Os` base, newlib-nano, `-ffast-math -fno-finite-math-only`; selective `HS_O3` column draw in `core/render/scan/shader.h` and pullback shader/pipeline paths |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | JewelMelt 288×144, single-entry playlist, source `6ab24361ce9699afe9b9365c0c7b6f360abde12d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 64 s capture; default one-hour epoch, no dwell compression |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261112-ashcloud-itcm HS_DEVICE_WAIT=1800 HS_TEENSY_PORT=COM4 HS_PROFILE_OUT=build/prof/jewelmelt_ship41b.log bash tools/profile_one.sh JewelMelt profile 64 32` |

Image size (measured single-entry Profile image):

```text
teensy_size:   FLASH: code:105888, data:151776, headers:8576   free for files:1765376
teensy_size:    RAM1: variables:314144, code:29992, padding:2776   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Shipping Phantasm image (41 entries):

```text
teensy_size:   FLASH: code:1183144, data:792876, headers:8492   free for files:47104
teensy_size:    RAM1: variables:313760, code:195656, padding:952   free for local variables:13920
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

The archived Phantasm build passes teensy-gate and hot-layout: RAM1 code is 195,656 B of 196,608 B, with 952 B remaining in six ITCM banks. RAM1 variables are 313,760 B; free local-variable space is 13,920 B against the 12,288 B stack floor.

Profile ELF SHA-256: `666b8da21d9aa17e8e24f60e003a9fa37882bfce5c9389fdb52688e0df05cdc7`. Phantasm ELF SHA-256: `050da474f94d77520d8f320508fee88272d55c44764393466064e3158204a162`. Both archived ELFs in `build/prof/artifacts/jewelmelt_ship41b_6ab24361ce96_666b8da21d9a` were rehashed for this report. Adjacent `jewelmelt_ship41b.provenance`, `jewelmelt_ship41b_envdump.txt` and `jewelmelt_ship41b_phantasm_envdump.txt` retain source, flags and compiler/ABI: Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 20251203.

Exactness cross-check: window frames 321–352, root 1,197,254,539 cyc ÷ 600 MHz = 1995424.232 us versus measured wall sum 1,995,428 us, within **1.888 ppm** (`python tools/parse_profile.py build/prof/jewelmelt_ship41b.log validate`; 5 ppm limit).

Hot-layout line in the measured Profile build (PASS):

```text
  ok   JewelMelt                worst 2 lines/set (0 sets over)  gated 20,478 B  overflow  5,216 B  full-set worst 2
```

Corresponding line in the 41-entry Phantasm build (PASS):

```text
  ok   JewelMelt                worst 2 lines/set (0 sets over)  gated 20,476 B  overflow  6,106 B  full-set worst 3
```

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer, setup included): peak frame render 79.887 ms, spilled 1/992 frames (0.10%). `fx_shader_draw` averages 37.403 ms/f; its worst window is 39.446 ms/f (frames 321–352).

**Steady render, frame 1 excluded:** peak **52.303 ms** (frame 317), **0/991 spilled** (0.00%). Setup frame 1: 79.887 ms. Complete windows cover frames 1–992; steady statistics use frames 2–992. Trailing rows without a complete counter window are excluded; subsequent transitions remain included. Steady mean render is 43.612 ms/frame.

A display window is 62.5 ms at 480 RPM; the effect renders one quadrant ≈ 10,368 px. All measured steady frames fit 16 fps, with 10.197 ms peak margin; steady render spans 38.492–52.303 ms. Setup exceeds one display window and is reported separately. `canvas_buffer_wait` is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the default shipping timeline runs through startup and a continuously evolving shader; no preset markers occur. The windows below show startup, the steady peak, the worst shader window and the lowest-cost steady window where distinct. Counter values are per-frame cycle-derived averages; percentages use the root frame, and nested costs overlap.

### Startup (window frames 1–32)

```text
frame                          62.273 ms  37.364 Mcyc 100.0%
  pov_preserve_half           134.378 us  80.627 kcyc   0.2%
    x0.97/frame; 138.713 us/call
  fx_shader_draw               38.098 ms  22.859 Mcyc  61.2%
    fx_prepare_column         729.208 us 437.525 kcyc   1.2%
      x144.00/frame; 5.064 us/call
  fx_palette_base_bake        846.412 us 507.847 kcyc   1.4%
    fx_palette_rgb_bake       720.826 us 432.495 kcyc   1.2%
      x1.00/frame; 720.826 us/call
  fx_prepare_frame              3.725 ms   2.235 Mcyc   6.0%
    x1.00/frame; 3724.909 us/call
  fx_advance                  982.669 us 589.602 kcyc   1.6%
    palette_cycle_step        952.684 us 571.611 kcyc   1.5%
      x1.00/frame; 952.684 us/call
  fx_timeline_step             85.445 us  51.267 kcyc   0.1%
    x1.00/frame; 85.445 us/call
  canvas_clear                 92.242 us  55.345 kcyc   0.1%
    x1.00/frame; 92.242 us/call
  canvas_buffer_wait           18.163 ms  10.898 Mcyc  29.2%
    x1.00/frame; 18162.840 us/call
```

Wall min/avg/max = 43.108/62.273/79.890 ms. Window mean render is 44.110 ms. This window includes setup frame 1; its maximum is excluded from the steady ranking. Later windows retain the changing shader workload without that setup draw.

### Peak-render regime (window frames 289–320)

```text
frame                          62.429 ms  37.458 Mcyc 100.0%
  pov_preserve_half           137.986 us  82.792 kcyc   0.2%
    x1.00/frame; 137.986 us/call
  fx_shader_draw               38.559 ms  23.136 Mcyc  61.8%
    fx_prepare_column         737.024 us 442.215 kcyc   1.2%
      x144.00/frame; 5.118 us/call
  fx_palette_base_bake          2.335 ms   1.401 Mcyc   3.7%
    fx_palette_rgb_bake         2.210 ms   1.326 Mcyc   3.5%
      x1.00/frame; 2209.841 us/call
  fx_prepare_frame              3.780 ms   2.268 Mcyc   6.1%
    x1.00/frame; 3779.604 us/call
  fx_advance                  895.238 us 537.143 kcyc   1.4%
    palette_cycle_step        862.723 us 517.634 kcyc   1.4%
      x1.00/frame; 862.723 us/call
  fx_timeline_step             85.059 us  51.035 kcyc   0.1%
    x1.00/frame; 85.059 us/call
  canvas_clear                 87.549 us  52.529 kcyc   0.1%
    x1.00/frame; 87.549 us/call
  canvas_buffer_wait           16.401 ms   9.840 Mcyc  26.3%
    x1.00/frame; 16400.593 us/call
```

Wall min/avg/max = 53.055/62.429/72.005 ms. Window mean render is 46.029 ms. This window contains the worst steady frame; the counter means do not equal that single-frame peak. The shader dominates live work, while the buffer wait fills the remaining display interval.

### Worst dominant-scope regime (window frames 321–352)

```text
frame                          62.357 ms  37.414 Mcyc 100.0%
  pov_preserve_half           138.692 us  83.215 kcyc   0.2%
    x1.00/frame; 138.692 us/call
  fx_shader_draw               39.446 ms  23.668 Mcyc  63.3%
    fx_prepare_column         754.317 us 452.590 kcyc   1.2%
      x144.00/frame; 5.238 us/call
  fx_palette_base_bake          2.605 ms   1.563 Mcyc   4.2%
    fx_palette_rgb_bake         2.480 ms   1.488 Mcyc   4.0%
      x1.00/frame; 2480.087 us/call
  fx_prepare_frame              3.723 ms   2.234 Mcyc   6.0%
    x1.00/frame; 3723.207 us/call
  fx_advance                  913.177 us 547.906 kcyc   1.5%
    palette_cycle_step        881.521 us 528.913 kcyc   1.4%
      x1.00/frame; 881.521 us/call
  fx_timeline_step             86.877 us  52.126 kcyc   0.1%
    x1.00/frame; 86.877 us/call
  canvas_clear                 87.304 us  52.383 kcyc   0.1%
    x1.00/frame; 87.304 us/call
  canvas_buffer_wait           15.211 ms   9.126 Mcyc  24.4%
    x1.00/frame; 15210.632 us/call
```

Wall min/avg/max = 53.007/62.357/71.836 ms. Window mean render is 47.146 ms. This window has the largest mean shader cost, although its single-frame render peak is lower than the peak-render window. The shader and preparation scopes explain the reduced idle interval relative to the lowest-cost regime.

### Lowest-cost steady regime (window frames 705–736)

```text
frame                          62.376 ms  37.426 Mcyc 100.0%
  pov_preserve_half           139.079 us  83.447 kcyc   0.2%
    x1.00/frame; 139.079 us/call
  fx_shader_draw               35.197 ms  21.118 Mcyc  56.4%
    fx_prepare_column         735.015 us 441.009 kcyc   1.2%
      x144.00/frame; 5.104 us/call
  fx_palette_base_bake        732.508 us 439.505 kcyc   1.2%
    fx_palette_rgb_bake       610.141 us 366.085 kcyc   1.0%
      x1.00/frame; 610.141 us/call
  fx_prepare_frame              3.531 ms   2.119 Mcyc   5.7%
    x1.00/frame; 3531.418 us/call
  fx_advance                  909.339 us 545.603 kcyc   1.5%
    palette_cycle_step        884.980 us 530.988 kcyc   1.4%
      x1.00/frame; 884.980 us/call
  fx_timeline_step             82.301 us  49.381 kcyc   0.1%
    x1.00/frame; 82.301 us/call
  canvas_clear                 85.848 us  51.509 kcyc   0.1%
    x1.00/frame; 85.848 us/call
  canvas_buffer_wait           21.571 ms  12.943 Mcyc  34.6%
    x1.00/frame; 21571.090 us/call
```

Wall min/avg/max = 57.579/62.376/67.101 ms. Window mean render is 40.805 ms. This is the lowest mean-render window after startup. Its larger buffer wait reflects the lower live workload at the same 16 fps cadence.

### Per-pixel figures

Per-pixel blend counters are disabled in this shipping timing image. Blended px/frame, quadrant coverage, cycles/blend and scan cycles/blended pixel are unavailable; geometric quadrant size is not substituted for measured blend counts.

## Column-ISR / DMA marshaling cost

Last complete window (frames 961–992); separately measured ISR interval = 2,001,420 us. Rates divide counts by 32 frames; CPU shares use that ISR interval, which includes report output.

```text
isr_wake        1152.88/frame
  min/avg/max 0.328/1.710/20.780 us  cpu 3.15%
  isr_pack        144.09/frame
    min/avg/max 7.576/8.313/14.420 us  cpu 1.91%
  isr_dma_submit  144.09/frame
    min/avg/max 0.746/0.918/9.195 us  cpu 0.21%
```

- CPU submission averages 0.918 us/call versus 8.313 us/call for packing (0.21% versus 1.91% CPU).
- Submission measures CPU work; asynchronous wire-transfer time is not measured.
- Inclusive `isr_wake` CPU share is 3.15%, leaving approximately 60.531 ms of foreground CPU per 62.5 ms window. Pack and submit run inside wake; use wake alone. DMA-completion and other interrupts are not measured. Render scopes already absorb interrupts, so no second ISR subtraction applies to their 62.5 ms budget. Every measured steady regime requires 1.000× speedup (already fits); setup alone would require 1.278× to fit one window.

## Summary ranking

1. `fx_shader_draw` — 61.76% of the frame, 38.559 ms: column shader draw.
2. `fx_prepare_frame` — 6.05% of the frame, 3.780 ms: frame preparation.
3. `fx_palette_base_bake` — 3.74% of the frame, 2.335 ms: palette base bake.
4. `fx_palette_rgb_bake` — 3.54% of the frame, 2.210 ms: RGB palette bake nested in base bake.
5. `fx_advance` — 1.43% of the frame, 0.895 ms: animation and palette advancement.

Ranked costs use the peak-render window, frames 289–320; nested scopes are inclusive. Native/WASM results are not a matched hardware timing baseline; this capture supplies no paired native/WASM measurement.

README cells: peak 🟢 52.303, spilled 🟢 0/991 (0.00%).

README setup frame 1: 79.887 ms. Captured 2026-11-13 07:02 America/Los_Angeles.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` first-caller parenting can hide its subtree when that parent has no calls; calls approximate blended pixels when enabled. Per-pixel scope overhead applies when enabled; per-pixel counters are disabled here.
- Shipping selective `HS_O3` column-draw and pullback shader/pipeline paths are active; this is the selective-O3 shipping configuration.
- No dwell-compression, ordered-cycle or transition-speed overrides were used; the default one-hour epoch exceeds this 64 s capture. These knobs control coverage rather than per-frame cost.
- Provenance records source `6ab24361ce9699afe9b9365c0c7b6f360abde12d`; working-tree cleanliness is not recorded in the run artifacts.
- This segment-0 master capture bounds the measured default-control interval; other segment roles and exhaustive parameter combinations are not measured.

Parser validation:

```text
effect=JewelMelt  windows=31
  [PASS] captured >=3 windows (31)
  [PASS] every window names one effect and resolution (JewelMelt 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 31 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 31 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.9 ppm, frames 321-352)

VALID: JewelMelt
```

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=JewelMelt`, `HS_PROFILE_WINDOW=32`; `just profile JewelMelt` = build + flash + capture through the shared-device lock. The Setup command reproduces this capture's duration, port and output path.
