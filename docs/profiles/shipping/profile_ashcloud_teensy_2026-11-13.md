# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-11-13, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `build/prof/ashcloud_ship41b.log` in `C:/work/temp/Holosphere-20261112-ashcloud-itcm`. Captured 2026-11-13 07:04 America/Los_Angeles on COM4 (serial 6362390). Replaces `shipping/profile_ashcloud_teensy_2026-11-05.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4 serial 6362390, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: shipping `-Os` base, newlib-nano, `-ffast-math -fno-finite-math-only`; selective `HS_O3` column draw in `core/render/scan/shader.h` and pullback shader/pipeline paths |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, source `6ab24361ce9699afe9b9365c0c7b6f360abde12d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 64 s capture; default one-hour epoch, no dwell compression |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261112-ashcloud-itcm HS_DEVICE_WAIT=1800 HS_TEENSY_PORT=COM4 HS_PROFILE_OUT=build/prof/ashcloud_ship41b.log bash tools/profile_one.sh AshCloud profile 64 32` |

Image size (measured single-entry Profile image):

```text
teensy_size:   FLASH: code:125424, data:182120, headers:8872   free for files:1715200
teensy_size:    RAM1: variables:314144, code:39656, padding:25880   free for local variables:144608
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Shipping Phantasm image (41 entries):

```text
teensy_size:   FLASH: code:1183144, data:792876, headers:8492   free for files:47104
teensy_size:    RAM1: variables:313760, code:195656, padding:952   free for local variables:13920
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

The archived Phantasm build passes teensy-gate and hot-layout: RAM1 code is 195,656 B of 196,608 B, with 952 B remaining in six ITCM banks. RAM1 variables are 313,760 B; free local-variable space is 13,920 B against the 12,288 B stack floor.

Profile ELF SHA-256: `5846bf08fef645a901d3df13f95b5f023538676d80ded561052fbacc6be78d0b`. Phantasm ELF SHA-256: `050da474f94d77520d8f320508fee88272d55c44764393466064e3158204a162`. Both archived ELFs in `build/prof/artifacts/ashcloud_ship41b_6ab24361ce96_5846bf08fef6` were rehashed for this report. Adjacent `ashcloud_ship41b.provenance`, `ashcloud_ship41b_envdump.txt` and `ashcloud_ship41b_phantasm_envdump.txt` retain source, flags and compiler/ABI: Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 20251203.

Exactness cross-check: window frames 897–928, root 1,198,281,380 cyc ÷ 600 MHz = 1997135.633 us versus measured wall sum 1,997,133 us, within **1.319 ppm** (`python tools/parse_profile.py build/prof/ashcloud_ship41b.log validate`; 5 ppm limit).

Hot-layout line in the measured Profile build (PASS):

```text
  ok   AshCloud                 worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 13,718 B  full-set worst 3
```

Corresponding line in the 41-entry Phantasm build (PASS):

```text
  ok   AshCloud                 worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 14,014 B  full-set worst 4
```

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer, setup included): peak frame render 85.409 ms, spilled 1/992 frames (0.10%). `fx_shader_draw` averages 45.745 ms/f; its worst window is 48.619 ms/f (frames 897–928).

**Steady render, frame 1 excluded:** peak **55.703 ms** (frame 920), **0/991 spilled** (0.00%). Setup frame 1: 85.409 ms. Complete windows cover frames 1–992; steady statistics use frames 2–992. Trailing rows without a complete counter window are excluded; subsequent transitions remain included. Steady mean render is 49.295 ms/frame.

A display window is 62.5 ms at 480 RPM; the effect renders one quadrant ≈ 10,368 px. All measured steady frames fit 16 fps, with 6.797 ms peak margin; steady render spans 38.221–55.703 ms. Setup exceeds one display window and is reported separately. `canvas_buffer_wait` is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the default shipping timeline runs through startup and a continuously evolving shader; no preset markers occur. The windows below show startup, the steady peak, the worst shader window and the lowest-cost steady window where distinct. Counter values are per-frame cycle-derived averages; percentages use the root frame, and nested costs overlap.

### Startup (window frames 1–32)

```text
frame                          62.498 ms  37.499 Mcyc 100.0%
  pov_preserve_half           135.980 us  81.588 kcyc   0.2%
    x0.97/frame; 140.366 us/call
  fx_shader_draw               41.246 ms  24.748 Mcyc  66.0%
    fx_prepare_column           2.106 ms   1.264 Mcyc   3.4%
      x144.00/frame; 14.627 us/call
  fx_palette_base_bake        570.278 us 342.167 kcyc   0.9%
    fx_palette_rgb_bake       451.133 us 270.680 kcyc   0.7%
      x1.00/frame; 451.133 us/call
  fx_prepare_frame              1.126 ms 675.418 kcyc   1.8%
    x1.00/frame; 1125.697 us/call
  fx_advance                    1.065 ms 638.967 kcyc   1.7%
    palette_cycle_step          1.014 ms 608.434 kcyc   1.6%
      x1.00/frame; 1014.056 us/call
  fx_timeline_step            125.116 us  75.070 kcyc   0.2%
    x1.00/frame; 125.116 us/call
  canvas_clear                 97.995 us  58.797 kcyc   0.2%
    x1.00/frame; 97.995 us/call
  canvas_buffer_wait           17.959 ms  10.775 Mcyc  28.7%
    x1.00/frame; 17958.680 us/call
```

Wall min/avg/max = 38.296/62.498/85.412 ms. Window mean render is 44.539 ms. This window includes setup frame 1; its maximum is excluded from the steady ranking. Later windows retain the changing shader workload without that setup draw.

### Peak-render regime / worst dominant-scope window (window frames 897–928)

```text
frame                          62.410 ms  37.446 Mcyc 100.0%
  pov_preserve_half           135.207 us  81.124 kcyc   0.2%
    x1.00/frame; 135.207 us/call
  fx_shader_draw               48.619 ms  29.171 Mcyc  77.9%
    fx_prepare_column           2.038 ms   1.223 Mcyc   3.3%
      x144.00/frame; 14.150 us/call
  fx_palette_base_bake        945.954 us 567.572 kcyc   1.5%
    fx_palette_rgb_bake       801.980 us 481.188 kcyc   1.3%
      x1.00/frame; 801.980 us/call
  fx_prepare_frame              1.279 ms 767.280 kcyc   2.0%
    x1.00/frame; 1278.799 us/call
  fx_advance                    1.028 ms 616.987 kcyc   1.6%
    palette_cycle_step        962.841 us 577.705 kcyc   1.5%
      x1.00/frame; 962.841 us/call
  fx_timeline_step            125.861 us  75.517 kcyc   0.2%
    x1.00/frame; 125.861 us/call
  canvas_clear                 96.486 us  57.891 kcyc   0.2%
    x1.00/frame; 96.486 us/call
  canvas_buffer_wait           10.015 ms   6.009 Mcyc  16.0%
    x1.00/frame; 10014.796 us/call
```

Wall min/avg/max = 57.698/62.410/67.204 ms. Window mean render is 52.396 ms. This window contains the worst steady frame; the counter means do not equal that single-frame peak. The shader dominates live work, while the buffer wait fills the remaining display interval.

### Lowest-cost steady regime (window frames 65–96)

```text
frame                          62.317 ms  37.390 Mcyc 100.0%
  pov_preserve_half           136.772 us  82.063 kcyc   0.2%
    x1.00/frame; 136.772 us/call
  fx_shader_draw               42.585 ms  25.551 Mcyc  68.3%
    fx_prepare_column           2.126 ms   1.276 Mcyc   3.4%
      x144.00/frame; 14.767 us/call
  fx_palette_base_bake        624.009 us 374.406 kcyc   1.0%
    fx_palette_rgb_bake       497.374 us 298.424 kcyc   0.8%
      x1.00/frame; 497.374 us/call
  fx_prepare_frame              1.279 ms 767.412 kcyc   2.1%
    x1.00/frame; 1279.019 us/call
  fx_advance                    1.035 ms 621.133 kcyc   1.7%
    palette_cycle_step        984.094 us 590.457 kcyc   1.6%
      x1.00/frame; 984.094 us/call
  fx_timeline_step            130.390 us  78.234 kcyc   0.2%
    x1.00/frame; 130.390 us/call
  canvas_clear                 94.342 us  56.605 kcyc   0.2%
    x1.00/frame; 94.342 us/call
  canvas_buffer_wait           16.261 ms   9.756 Mcyc  26.1%
    x1.00/frame; 16260.702 us/call
```

Wall min/avg/max = 58.863/62.317/66.290 ms. Window mean render is 46.056 ms. This is the lowest mean-render window after startup. Its larger buffer wait reflects the lower live workload at the same 16 fps cadence.

### Per-pixel figures

Per-pixel blend counters are disabled in this shipping timing image. Blended px/frame, quadrant coverage, cycles/blend and scan cycles/blended pixel are unavailable; geometric quadrant size is not substituted for measured blend counts.

## Column-ISR / DMA marshaling cost

Last complete window (frames 961–992); separately measured ISR interval = 1,999,700 us. Rates divide counts by 32 frames; CPU shares use that ISR interval, which includes report output.

```text
isr_wake        1151.91/frame
  min/avg/max 0.285/1.768/23.121 us  cpu 3.26%
  isr_pack        143.97/frame
    min/avg/max 8.053/8.708/11.501 us  cpu 2.00%
  isr_dma_submit  143.97/frame
    min/avg/max 0.596/0.918/7.113 us  cpu 0.21%
```

- CPU submission averages 0.918 us/call versus 8.708 us/call for packing (0.21% versus 2.00% CPU).
- Submission measures CPU work; asynchronous wire-transfer time is not measured.
- Inclusive `isr_wake` CPU share is 3.26%, leaving approximately 60.463 ms of foreground CPU per 62.5 ms window. Pack and submit run inside wake; use wake alone. DMA-completion and other interrupts are not measured. Render scopes already absorb interrupts, so no second ISR subtraction applies to their 62.5 ms budget. Every measured steady regime requires 1.000× speedup (already fits); setup alone would require 1.367× to fit one window.

## Summary ranking

1. `fx_shader_draw` — 77.90% of the frame, 48.619 ms: column shader draw.
2. `fx_prepare_column` — 3.26% of the frame, 2.038 ms: column preparation nested in shader draw.
3. `fx_prepare_frame` — 2.05% of the frame, 1.279 ms: frame preparation.
4. `fx_advance` — 1.65% of the frame, 1.028 ms: animation and palette advancement.
5. `palette_cycle_step` — 1.54% of the frame, 0.963 ms: palette cycling nested in advancement.

Ranked costs use the peak-render window, frames 897–928; nested scopes are inclusive. Native/WASM results are not a matched hardware timing baseline; this capture supplies no paired native/WASM measurement.

README cells: peak 🟢 55.703, spilled 🟢 0/991 (0.00%).

README setup frame 1: 85.409 ms. Captured 2026-11-13 07:04 America/Los_Angeles.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` first-caller parenting can hide its subtree when that parent has no calls; calls approximate blended pixels when enabled. Per-pixel scope overhead applies when enabled; per-pixel counters are disabled here.
- Shipping selective `HS_O3` column-draw and pullback shader/pipeline paths are active; this is the selective-O3 shipping configuration.
- No dwell-compression, ordered-cycle or transition-speed overrides were used; the default one-hour epoch exceeds this 64 s capture. These knobs control coverage rather than per-frame cost.
- Provenance records source `6ab24361ce9699afe9b9365c0c7b6f360abde12d`; working-tree cleanliness is not recorded in the run artifacts.
- Source carries AshCloud fallback commit `b77c0835c`, which outlines cached shading fallbacks into flash.
- This segment-0 master capture bounds the measured default-control interval; other segment roles and exhaustive parameter combinations are not measured.

Parser validation:

```text
effect=AshCloud  windows=31
  [PASS] captured >=3 windows (31)
  [PASS] every window names one effect and resolution (AshCloud 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 31 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 31 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.3 ppm, frames 897-928)

VALID: AshCloud
```

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`; `just profile AshCloud` = build + flash + capture through the shared-device lock. The Setup command reproduces this capture's duration, port and output path.
