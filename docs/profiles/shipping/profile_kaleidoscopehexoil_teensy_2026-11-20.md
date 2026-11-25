# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 105-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/kaleidoscopehexoil_ship.log`. Supersedes `profile_kaleidoscopehexoil_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | KaleidoscopeHexOil at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 105 s slot, 111 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh KaleidoscopeHexOil profile 111 16` from the captured revision |

Captured 2026-11-20 15:47, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `68353b994b7e3ddb4417127c97923f9c5b38145c390d46fe80f15c0ebf474302`.

Profile image size `FLASH: code:103488, data:156568, headers:8232 free for files:1763328 / RAM1: variables:314144, code:26200, padding:6568 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `540bb01f2c65bcbe5320f97eaca8ca9ea0926c7aae356ef4733a272aa33290a4`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 849–864. The root counter measured 725395.425 µs, compared with 725,393 µs of draw time (3.34 ppm).

## Frame cadence

**Peak confirmed live render 49.700 ms** at frame 882 (kaleidoscope-hex-oil → kaleidoscope-hex-oil-2). **0/1,668 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 71.154 ms. Live render averaged 40.214 ms.

Among complete live readouts, the largest window mean was 45.336 ms/frame at frames 849–864. The dominant scope averaged 37.086 ms/frame across those readouts.

The shipping ELF assigns 840 revolutions (105 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 105 s. Frames 2–1669 are confirmed live. The raw file contains 1,752 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 10 frame(s), 1670–1679, with uncertain membership. Their largest render was 38.570 ms and 0 exceeded the budget. The live denominator lies between 1,668 and 1,678. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 49.700 (2), spilled 🟢 0/1668 (0.00%).

README setup frame 1: 71.154 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    39.774 ms 23.86 Mcyc 100.0%
  pov_preserve_half      130.30 us  78.2 kcyc   0.3% x0.94 139 us/call
  fx_shader_draw         37.156 ms 22.29 Mcyc  93.4% x1 37156 us/call
    fx_prepare_column    772.27 us 463.4 kcyc   1.9% x144 5.36 us/call
  fx_palette_base_bake   637.55 us 382.5 kcyc   1.6% x1 638 us/call
    fx_palette_rgb_bake  508.83 us 305.3 kcyc   1.3% x1 509 us/call
  fx_prepare_frame       549.08 us 329.5 kcyc   1.4% x1 549 us/call
  fx_advance             954.58 us 572.7 kcyc   2.4% x1 955 us/call
    palette_cycle_step   911.36 us 546.8 kcyc   2.3% x1 911 us/call
  fx_timeline_step        78.85 us  47.3 kcyc   0.2% x1 78.8 us/call
  canvas_clear           101.99 us  61.2 kcyc   0.3% x1 102 us/call
  canvas_buffer_wait       2.88 us   1.7 kcyc   0.0% x1 2.88 us/call
```

Wall min/mean/max was 36.255/39.774/71.156 ms. Render averaged 39.772 ms per frame, including 37.156 ms in `fx_shader_draw`.

This regime leaves 22.728 ms of mean render headroom per display window.

### Window containing the live peak (frames 881–896)

```
frame                    45.321 ms 27.19 Mcyc 100.0%
  pov_preserve_half      137.95 us  82.8 kcyc   0.3% x1 138 us/call
  fx_shader_draw         41.014 ms 24.61 Mcyc  90.5% x1 41014 us/call
    fx_prepare_column     1.951 ms  1.17 Mcyc   4.3% x144 13.5 us/call
  fx_palette_base_bake    2.139 ms  1.28 Mcyc   4.7% x1 2139 us/call
    fx_palette_rgb_bake   2.000 ms  1.20 Mcyc   4.4% x1 2000 us/call
  fx_prepare_frame       602.40 us 361.4 kcyc   1.3% x1 602 us/call
  fx_advance             983.87 us 590.3 kcyc   2.2% x1 984 us/call
    palette_cycle_step   946.58 us 567.9 kcyc   2.1% x1 947 us/call
  fx_timeline_step       173.83 us 104.3 kcyc   0.4% x1 174 us/call
  canvas_clear            95.64 us  57.4 kcyc   0.2% x1 95.6 us/call
  canvas_buffer_wait       1.40 us    841 cyc   0.0% x1 1.40 us/call
```

Wall min/mean/max was 40.590/45.320/49.701 ms. Render averaged 45.320 ms per frame, including 41.014 ms in `fx_shader_draw`.

This regime leaves 17.180 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 353–368)

```
frame                    41.735 ms 25.04 Mcyc 100.0%
  pov_preserve_half      139.64 us  83.8 kcyc   0.3% x1 140 us/call
  fx_shader_draw         38.300 ms 22.98 Mcyc  91.8% x1 38300 us/call
    fx_prepare_column    777.05 us 466.2 kcyc   1.9% x144 5.40 us/call
  fx_palette_base_bake    1.332 ms 798.9 kcyc   3.2% x1 1332 us/call
    fx_palette_rgb_bake   1.199 ms 719.7 kcyc   2.9% x1 1199 us/call
  fx_prepare_frame       690.59 us 414.4 kcyc   1.7% x1 691 us/call
  fx_advance             954.24 us 572.5 kcyc   2.3% x1 954 us/call
    palette_cycle_step   913.04 us 547.8 kcyc   2.2% x1 913 us/call
  fx_timeline_step        70.54 us  42.3 kcyc   0.2% x1 70.5 us/call
  canvas_clear            94.56 us  56.7 kcyc   0.2% x1 94.6 us/call
  canvas_buffer_wait       2.47 us   1.5 kcyc   0.0% x1 2.47 us/call
```

Wall min/mean/max was 40.300/41.735/43.360 ms. Render averaged 41.733 ms per frame, including 38.300 ms in `fx_shader_draw`.

This regime leaves 20.767 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 849–864)

```
frame                    45.337 ms 27.20 Mcyc 100.0%
  pov_preserve_half      137.20 us  82.3 kcyc   0.3% x1 137 us/call
  fx_shader_draw         41.115 ms 24.67 Mcyc  90.7% x1 41115 us/call
    fx_prepare_column     1.930 ms  1.16 Mcyc   4.3% x144 13.4 us/call
  fx_palette_base_bake    2.043 ms  1.23 Mcyc   4.5% x1 2043 us/call
    fx_palette_rgb_bake   1.907 ms  1.14 Mcyc   4.2% x1 1907 us/call
  fx_prepare_frame       612.62 us 367.6 kcyc   1.4% x1 613 us/call
  fx_advance             988.96 us 593.4 kcyc   2.2% x1 989 us/call
    palette_cycle_step   949.75 us 569.8 kcyc   2.1% x1 950 us/call
  fx_timeline_step       172.57 us 103.5 kcyc   0.4% x1 173 us/call
  canvas_clear            95.76 us  57.5 kcyc   0.2% x1 95.8 us/call
  canvas_buffer_wait       1.54 us    924 cyc   0.0% x1 1.54 us/call
```

Wall min/mean/max was 41.154/45.337/49.656 ms. Render averaged 45.336 ms per frame, including 41.115 ms in `fx_shader_draw`.

This regime leaves 17.164 ms of mean render headroom per display window.

### Final complete in-show window (frames 1649–1664)

```
frame                    38.834 ms 23.30 Mcyc 100.0%
  pov_preserve_half      139.28 us  83.6 kcyc   0.4% x1 139 us/call
  fx_shader_draw         35.834 ms 21.50 Mcyc  92.3% x1 35834 us/call
    fx_prepare_column    773.02 us 463.8 kcyc   2.0% x144 5.37 us/call
  fx_palette_base_bake   999.43 us 599.7 kcyc   2.6% x1 999 us/call
    fx_palette_rgb_bake  866.78 us 520.1 kcyc   2.2% x1 867 us/call
  fx_prepare_frame       583.16 us 349.9 kcyc   1.5% x1 583 us/call
  fx_advance             964.48 us 578.7 kcyc   2.5% x1 964 us/call
    palette_cycle_step   921.96 us 553.2 kcyc   2.4% x1 922 us/call
  fx_timeline_step        71.64 us  43.0 kcyc   0.2% x1 71.6 us/call
  canvas_clear            93.57 us  56.1 kcyc   0.2% x1 93.6 us/call
  canvas_buffer_wait       2.95 us   1.8 kcyc   0.0% x1 2.95 us/call
```

Wall min/mean/max was 38.331/38.833/39.603 ms. Render averaged 38.831 ms per frame, including 35.834 ms in `fx_shader_draw`.

This regime leaves 23.669 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| kaleidoscope-hex-oil | held preset | 36 | 353–368 | — | — | 38.300 | 41.733 | 47.737 | 16 | 0/598 |
| kaleidoscope-hex-oil-2 | held preset | 36 | 1297–1312 | — | — | 37.340 | 39.981 | 49.700 | 16 | 0/1070 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 881–896, measured over the separate 999,440 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.4/frame  0.581/1.783/19.453 us  3.28%
  isr_pack           143.9/frame  6.853/7.690/10.078 us  1.77%
  isr_dma_submit     143.9/frame  0.670/0.920/6.591 us  0.21%
```

Inclusive wake time leaves 60.450 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 41.014 ms/frame (90.5% of draw time).
2. `fx_palette_base_bake` took 2.139 ms/frame (4.7% of draw time).
3. `fx_palette_rgb_bake` took 2.000 ms/frame (4.4% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeHexOil` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
