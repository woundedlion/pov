# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 240-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/kaleidoscopesmooth_ship.log`. Supersedes `profile_kaleidoscopesmooth_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | KaleidoscopeSmooth at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 240 s slot, 246 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeSmooth profile 246 16` from the captured revision |

Captured 2026-11-20 15:10, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `6522c8033615a46773359c1c9364fcfd74ef90d615b66bfdf872d929e21c2af6`.

Profile image size `FLASH: code:102368, data:153008, headers:8816 free for files:1767424 / RAM1: variables:314144, code:27464, padding:5304 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `9202624cc128dfe99bd14537f61c2d9fb2f2d2bab0294454c53a529d44232427`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 833–848. The root counter measured 611882.125 µs, compared with 611,880 µs of draw time (3.47 ppm).

## Frame cadence

**Peak confirmed live render 41.435 ms** at frame 307 (coupled-grid). **0/3,825 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 65.479 ms. Live render averaged 32.508 ms.

Among complete live readouts, the largest window mean was 38.241 ms/frame at frames 833–848. The dominant scope averaged 29.104 ms/frame across those readouts.

The shipping ELF assigns 1,920 revolutions (240 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 240 s. Frames 2–3826 are confirmed live. The raw file contains 3,912 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 3827–3840, with uncertain membership. Their largest render was 36.935 ms and 0 exceeded the budget. The live denominator lies between 3,825 and 3,839. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 41.435 (4), spilled 🟢 0/3825 (0.00%).

README setup frame 1: 65.479 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    34.538 ms 20.72 Mcyc 100.0%
  pov_preserve_half      131.75 us  79.1 kcyc   0.4% x0.94 141 us/call
  fx_shader_draw         31.311 ms 18.79 Mcyc  90.7% x1 31311 us/call
    fx_prepare_column    627.96 us 376.8 kcyc   1.8% x144 4.36 us/call
  fx_palette_base_bake   804.47 us 482.7 kcyc   2.3% x1 804 us/call
    fx_palette_rgb_bake  667.41 us 400.4 kcyc   1.9% x1 667 us/call
  fx_prepare_frame       894.57 us 536.7 kcyc   2.6% x1 895 us/call
  fx_advance              1.009 ms 605.7 kcyc   2.9% x1 1009 us/call
    palette_cycle_step   956.93 us 574.2 kcyc   2.8% x1 957 us/call
  fx_timeline_step        86.23 us  51.7 kcyc   0.2% x1 86.2 us/call
  canvas_clear           103.84 us  62.3 kcyc   0.3% x1 104 us/call
  canvas_buffer_wait       1.65 us    990 cyc   0.0% x1 1.65 us/call
```

Wall min/mean/max was 30.095/34.537/65.479 ms. Render averaged 34.536 ms per frame, including 31.311 ms in `fx_shader_draw`.

This regime leaves 27.964 ms of mean render headroom per display window.

### Window containing the live peak (frames 305–320)

```
frame                    36.326 ms 21.80 Mcyc 100.0%
  pov_preserve_half      139.94 us  84.0 kcyc   0.4% x1 140 us/call
  fx_shader_draw         31.434 ms 18.86 Mcyc  86.5% x1 31434 us/call
    fx_prepare_column    654.37 us 392.6 kcyc   1.8% x144 4.54 us/call
  fx_palette_base_bake    2.630 ms  1.58 Mcyc   7.2% x1 2630 us/call
    fx_palette_rgb_bake   2.495 ms  1.50 Mcyc   6.9% x1 2495 us/call
  fx_prepare_frame       808.44 us 485.1 kcyc   2.2% x1 808 us/call
  fx_advance             950.18 us 570.1 kcyc   2.6% x1 950 us/call
    palette_cycle_step   902.59 us 541.6 kcyc   2.5% x1 903 us/call
  fx_timeline_step        81.67 us  49.0 kcyc   0.2% x1 81.7 us/call
  canvas_clear            95.82 us  57.5 kcyc   0.3% x1 95.8 us/call
  canvas_buffer_wait       1.93 us   1.2 kcyc   0.0% x1 1.93 us/call
```

Wall min/mean/max was 31.990/36.326/41.436 ms. Render averaged 36.325 ms per frame, including 31.434 ms in `fx_shader_draw`.

This regime leaves 26.175 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 81–96)

```
frame                    36.149 ms 21.69 Mcyc 100.0%
  pov_preserve_half      137.53 us  82.5 kcyc   0.4% x1 138 us/call
  fx_shader_draw         32.846 ms 19.71 Mcyc  90.9% x1 32846 us/call
    fx_prepare_column    654.48 us 392.7 kcyc   1.8% x144 4.54 us/call
  fx_palette_base_bake    1.083 ms 649.8 kcyc   3.0% x1 1083 us/call
    fx_palette_rgb_bake  942.89 us 565.7 kcyc   2.6% x1 943 us/call
  fx_prepare_frame       758.18 us 454.9 kcyc   2.1% x1 758 us/call
  fx_advance             957.96 us 574.8 kcyc   2.7% x1 958 us/call
    palette_cycle_step   909.81 us 545.9 kcyc   2.5% x1 910 us/call
  fx_timeline_step        81.63 us  49.0 kcyc   0.2% x1 81.6 us/call
  canvas_clear            95.15 us  57.1 kcyc   0.3% x1 95.1 us/call
  canvas_buffer_wait       2.10 us   1.3 kcyc   0.0% x1 2.10 us/call
```

Wall min/mean/max was 32.775/36.148/39.814 ms. Render averaged 36.147 ms per frame, including 32.846 ms in `fx_shader_draw`.

This regime leaves 26.353 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 833–848)

```
frame                    38.243 ms 22.95 Mcyc 100.0%
  pov_preserve_half      139.85 us  83.9 kcyc   0.4% x1 140 us/call
  fx_shader_draw         35.062 ms 21.04 Mcyc  91.7% x1 35062 us/call
    fx_prepare_column    817.86 us 490.7 kcyc   2.1% x144 5.68 us/call
  fx_palette_base_bake    1.001 ms 600.7 kcyc   2.6% x1 1001 us/call
    fx_palette_rgb_bake  859.82 us 515.9 kcyc   2.2% x1 860 us/call
  fx_prepare_frame       600.81 us 360.5 kcyc   1.6% x1 601 us/call
  fx_advance             957.59 us 574.6 kcyc   2.5% x1 958 us/call
    palette_cycle_step   903.62 us 542.2 kcyc   2.4% x1 904 us/call
  fx_timeline_step       196.56 us 117.9 kcyc   0.5% x1 197 us/call
  canvas_clear            96.02 us  57.6 kcyc   0.3% x1 96.0 us/call
  canvas_buffer_wait       1.90 us   1.1 kcyc   0.0% x1 1.90 us/call
```

Wall min/mean/max was 35.470/38.242/40.275 ms. Render averaged 38.241 ms per frame, including 35.062 ms in `fx_shader_draw`.

This regime leaves 24.259 ms of mean render headroom per display window.

### Final complete in-show window (frames 3809–3824)

```
frame                    29.970 ms 17.98 Mcyc 100.0%
  pov_preserve_half      143.62 us  86.2 kcyc   0.5% x1 144 us/call
  fx_shader_draw         27.127 ms 16.28 Mcyc  90.5% x1 27127 us/call
    fx_prepare_column    689.73 us 413.8 kcyc   2.3% x144 4.79 us/call
  fx_palette_base_bake   902.24 us 541.3 kcyc   3.0% x1 902 us/call
    fx_palette_rgb_bake  762.27 us 457.4 kcyc   2.5% x1 762 us/call
  fx_prepare_frame       534.58 us 320.7 kcyc   1.8% x1 535 us/call
  fx_advance             895.48 us 537.3 kcyc   3.0% x1 895 us/call
    palette_cycle_step   846.08 us 507.6 kcyc   2.8% x1 846 us/call
  fx_timeline_step        81.12 us  48.7 kcyc   0.3% x1 81.1 us/call
  canvas_clear            95.27 us  57.2 kcyc   0.3% x1 95.3 us/call
  canvas_buffer_wait       2.10 us   1.3 kcyc   0.0% x1 2.10 us/call
```

Wall min/mean/max was 28.514/29.969/31.428 ms. Render averaged 29.968 ms per frame, including 27.127 ms in `fx_shader_draw`.

This regime leaves 32.532 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| coupled-grid | held preset | 36 | 81–96 | — | — | 32.846 | 36.147 | 41.435 | 16 | 0/598 |
| direct-grid | held preset | 36 | 1569–1584 | — | — | 31.391 | 34.988 | 41.061 | 16 | 0/1079 |
| double-map | held preset | 37 | 2641–2656 | — | — | 31.133 | 34.540 | 39.798 | 16 | 0/1079 |
| stretched-grid | held preset | 36 | 3521–3536 | — | — | 29.044 | 33.540 | 38.713 | 16 | 0/1069 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 305–320, measured over the separate 998,332 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1150.1/frame  0.500/1.743/24.911 us  3.21%
  isr_pack           143.8/frame  6.898/7.470/10.108 us  1.72%
  isr_dma_submit     143.8/frame  0.668/0.916/2.961 us  0.21%
```

Inclusive wake time leaves 60.494 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 31.434 ms/frame (86.5% of draw time).
2. `fx_palette_base_bake` took 2.630 ms/frame (7.2% of draw time).
3. `fx_palette_rgb_bake` took 2.495 ms/frame (6.9% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeSmooth` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
