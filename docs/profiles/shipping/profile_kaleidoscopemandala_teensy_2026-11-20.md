# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 105-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/kaleidoscopemandala_ship.log`. Supersedes `profile_kaleidoscopemandala_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | KaleidoscopeMandala at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 105 s slot, 111 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh KaleidoscopeMandala profile 111 16` from the captured revision |

Captured 2026-11-20 16:07, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `9692e4627728c07048d95a8603a949af8ed1afee19ad883d6058c255826245e3`.

Profile image size `FLASH: code:102792, data:153376, headers:9048 free for files:1766400 / RAM1: variables:314144, code:27400, padding:5368 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `9778dc8c47de2dcb614381158eb946f89580c2b22ff661dfc5db9a5df29da12d`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 849–864. The root counter measured 701900.637 µs, compared with 701,902 µs of draw time (1.94 ppm).

## Frame cadence

**Peak confirmed live render 47.556 ms** at frame 862 (wave-mirror → cup-hue). **0/1,668 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 67.029 ms. Live render averaged 38.560 ms.

Among complete live readouts, the largest window mean was 43.866 ms/frame at frames 849–864. The dominant scope averaged 34.183 ms/frame across those readouts.

The shipping ELF assigns 840 revolutions (105 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 105 s. Frames 2–1669 are confirmed live. The raw file contains 1,752 frames, including 74 confirmed after the boundary.

Boundary resolution leaves 9 frame(s), 1670–1678, with uncertain membership. Their largest render was 40.755 ms and 0 exceeded the budget. The live denominator lies between 1,668 and 1,677. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 47.556 (2), spilled 🟢 0/1668 (0.00%).

README setup frame 1: 67.029 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    36.624 ms 21.97 Mcyc 100.0%
  pov_preserve_half      131.80 us  79.1 kcyc   0.4% x0.94 141 us/call
  fx_shader_draw         33.418 ms 20.05 Mcyc  91.2% x1 33418 us/call
    fx_prepare_column    690.53 us 414.3 kcyc   1.9% x144 4.80 us/call
  fx_palette_base_bake   698.32 us 419.0 kcyc   1.9% x1 698 us/call
    fx_palette_rgb_bake  573.00 us 343.8 kcyc   1.6% x1 573 us/call
  fx_prepare_frame       877.52 us 526.5 kcyc   2.4% x1 878 us/call
  fx_advance              1.092 ms 655.1 kcyc   3.0% x1 1092 us/call
    palette_cycle_step    1.056 ms 633.7 kcyc   2.9% x1 1056 us/call
  fx_timeline_step        67.19 us  40.3 kcyc   0.2% x1 67.2 us/call
  canvas_clear           102.29 us  61.4 kcyc   0.3% x1 102 us/call
  canvas_buffer_wait       2.61 us   1.6 kcyc   0.0% x1 2.61 us/call
```

Wall min/mean/max was 32.698/36.624/67.031 ms. Render averaged 36.622 ms per frame, including 33.418 ms in `fx_shader_draw`.

This regime leaves 25.878 ms of mean render headroom per display window.

### Window containing the live peak (frames 849–864)

```
frame                    43.869 ms 26.32 Mcyc 100.0%
  pov_preserve_half      140.92 us  84.6 kcyc   0.3% x1 141 us/call
  fx_shader_draw         36.548 ms 21.93 Mcyc  83.3% x1 36548 us/call
    fx_prepare_column     1.083 ms 649.6 kcyc   2.5% x144 7.52 us/call
  fx_palette_base_bake    2.318 ms  1.39 Mcyc   5.3% x1 2318 us/call
    fx_palette_rgb_bake   2.195 ms  1.32 Mcyc   5.0% x1 2195 us/call
  fx_prepare_frame        3.213 ms  1.93 Mcyc   7.3% x1 3213 us/call
  fx_advance              1.138 ms 682.8 kcyc   2.6% x1 1138 us/call
    palette_cycle_step    1.095 ms 657.3 kcyc   2.5% x1 1095 us/call
  fx_timeline_step       195.03 us 117.0 kcyc   0.4% x1 195 us/call
  canvas_clear            95.35 us  57.2 kcyc   0.2% x1 95.3 us/call
  canvas_buffer_wait       2.69 us   1.6 kcyc   0.0% x1 2.69 us/call
```

Wall min/mean/max was 40.299/43.868/47.559 ms. Render averaged 43.866 ms per frame, including 36.548 ms in `fx_shader_draw`.

This regime leaves 18.634 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 545–560)

```
frame                    39.532 ms 23.72 Mcyc 100.0%
  pov_preserve_half      140.33 us  84.2 kcyc   0.4% x1 140 us/call
  fx_shader_draw         35.887 ms 21.53 Mcyc  90.8% x1 35887 us/call
    fx_prepare_column    707.31 us 424.4 kcyc   1.8% x144 4.91 us/call
  fx_palette_base_bake    1.193 ms 716.0 kcyc   3.0% x1 1193 us/call
    fx_palette_rgb_bake   1.068 ms 641.0 kcyc   2.7% x1 1068 us/call
  fx_prepare_frame       842.91 us 505.7 kcyc   2.1% x1 843 us/call
  fx_advance              1.079 ms 647.3 kcyc   2.7% x1 1079 us/call
    palette_cycle_step    1.043 ms 625.8 kcyc   2.6% x1 1043 us/call
  fx_timeline_step        62.31 us  37.4 kcyc   0.2% x1 62.3 us/call
  canvas_clear            96.27 us  57.8 kcyc   0.2% x1 96.3 us/call
  canvas_buffer_wait       2.59 us   1.6 kcyc   0.0% x1 2.59 us/call
```

Wall min/mean/max was 36.947/39.532/41.739 ms. Render averaged 39.529 ms per frame, including 35.887 ms in `fx_shader_draw`.

This regime leaves 22.971 ms of mean render headroom per display window.

### Final complete in-show window (frames 1649–1664)

```
frame                    39.369 ms 23.62 Mcyc 100.0%
  pov_preserve_half      140.38 us  84.2 kcyc   0.4% x1 140 us/call
  fx_shader_draw         35.789 ms 21.47 Mcyc  90.9% x1 35789 us/call
    fx_prepare_column    698.61 us 419.2 kcyc   1.8% x144 4.85 us/call
  fx_palette_base_bake    1.236 ms 741.9 kcyc   3.1% x1 1236 us/call
    fx_palette_rgb_bake   1.108 ms 664.8 kcyc   2.8% x1 1108 us/call
  fx_prepare_frame       755.51 us 453.3 kcyc   1.9% x1 756 us/call
  fx_advance              1.059 ms 635.4 kcyc   2.7% x1 1059 us/call
    palette_cycle_step    1.022 ms 613.4 kcyc   2.6% x1 1022 us/call
  fx_timeline_step        67.22 us  40.3 kcyc   0.2% x1 67.2 us/call
  canvas_clear            94.59 us  56.8 kcyc   0.2% x1 94.6 us/call
  canvas_buffer_wait       2.50 us   1.5 kcyc   0.0% x1 2.50 us/call
```

Wall min/mean/max was 37.360/39.369/41.171 ms. Render averaged 39.367 ms per frame, including 35.789 ms in `fx_shader_draw`.

This regime leaves 23.133 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| wave-mirror | held preset | 36 | 545–560 | — | — | 35.887 | 39.529 | 45.093 | 16 | 0/598 |
| cup-hue | held preset | 36 | 1649–1664 | — | — | 35.789 | 39.367 | 47.556 | 16 | 0/1070 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 849–864, measured over the separate 1,000,615 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.8/frame  0.578/1.781/23.700 us  3.28%
  isr_pack           144.1/frame  7.036/7.630/10.050 us  1.75%
  isr_dma_submit     144.1/frame  0.730/0.920/7.653 us  0.21%
```

Inclusive wake time leaves 60.450 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 36.548 ms/frame (83.3% of draw time).
2. `fx_prepare_frame` took 3.213 ms/frame (7.3% of draw time).
3. `fx_palette_base_bake` took 2.318 ms/frame (5.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeMandala` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
