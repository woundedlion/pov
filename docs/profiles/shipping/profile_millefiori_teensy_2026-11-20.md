# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/millefiori_ship.log`. Supersedes `profile_millefiori_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Millefiori at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Millefiori profile 66 32` from the captured revision |

Captured 2026-11-20 16:04, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `e0d4bcf1b6cf4318c830a712ed7e483e312cf4044e2a8bd55e6f1bb47babeca4`.

Profile image size `FLASH: code:108904, data:153144, headers:8288 free for files:1761280 / RAM1: variables:314112, code:28424, padding:4344 free for local variables:177408 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `50a25934427f69738fb81abcd89a8873e4ce6026147ce6d954ec77538a416b15`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 193–224. The root counter measured 1653589.277 µs, compared with 1,653,591 µs of draw time (1.04 ppm).

## Frame cadence

**Peak confirmed live render 53.861 ms** at frame 385 (Millefiori). **0/950 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 82.441 ms. Live render averaged 48.576 ms.

Among complete live readouts, the largest window mean was 51.673 ms/frame at frames 193–224. The dominant scope averaged 46.145 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–951 are confirmed live. The raw file contains 1,032 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 8 frame(s), 952–959, with uncertain membership. Their largest render was 53.074 ms and 0 exceeded the budget. The live denominator lies between 950 and 958. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 53.861, spilled 🟢 0/950 (0.00%).

README setup frame 1: 82.441 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    47.658 ms 28.59 Mcyc 100.0%
  pov_preserve_half      133.03 us  79.8 kcyc   0.3% x0.97 137 us/call
  fx_shader_draw         44.997 ms 27.00 Mcyc  94.4% x1 44997 us/call
    fx_prepare_column     2.474 ms  1.48 Mcyc   5.2% x144 17.2 us/call
  fx_palette_base_bake   116.37 us  69.8 kcyc   0.2% x1 116 us/call
  fx_prepare_frame       926.54 us 555.9 kcyc   1.9% x1 927 us/call
  fx_advance              1.099 ms 659.4 kcyc   2.3% x1 1099 us/call
    palette_cycle_step    1.044 ms 626.6 kcyc   2.2% x1 1044 us/call
  fx_timeline_step       107.37 us  64.4 kcyc   0.2% x1 107 us/call
  canvas_clear            90.77 us  54.5 kcyc   0.2% x1 90.8 us/call
  canvas_buffer_wait       1.98 us   1.2 kcyc   0.0% x1 1.98 us/call
```

Wall min/mean/max was 41.691/47.658/82.442 ms. Render averaged 47.656 ms per frame, including 44.997 ms in `fx_shader_draw`.

This regime leaves 14.844 ms of mean render headroom per display window.

### Window containing the live peak (frames 385–416)

```
frame                    48.767 ms 29.26 Mcyc 100.0%
  pov_preserve_half      137.52 us  82.5 kcyc   0.3% x1 138 us/call
  fx_shader_draw         46.370 ms 27.82 Mcyc  95.1% x1 46370 us/call
    fx_prepare_column     2.499 ms  1.50 Mcyc   5.1% x144 17.4 us/call
  fx_palette_base_bake   121.96 us  73.2 kcyc   0.3% x1 122 us/call
  fx_prepare_frame       827.22 us 496.3 kcyc   1.7% x1 827 us/call
  fx_advance             924.55 us 554.7 kcyc   1.9% x1 925 us/call
    palette_cycle_step   869.85 us 521.9 kcyc   1.8% x1 870 us/call
  fx_timeline_step       108.08 us  64.8 kcyc   0.2% x1 108 us/call
  canvas_clear            86.28 us  51.8 kcyc   0.2% x1 86.3 us/call
  canvas_buffer_wait       1.89 us   1.1 kcyc   0.0% x1 1.89 us/call
```

Wall min/mean/max was 45.305/48.767/53.862 ms. Render averaged 48.766 ms per frame, including 46.370 ms in `fx_shader_draw`.

This regime leaves 13.734 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    48.227 ms 28.94 Mcyc 100.0%
  pov_preserve_half      138.23 us  82.9 kcyc   0.3% x1 138 us/call
  fx_shader_draw         45.778 ms 27.47 Mcyc  94.9% x1 45778 us/call
    fx_prepare_column     2.489 ms  1.49 Mcyc   5.2% x144 17.3 us/call
  fx_palette_base_bake   118.89 us  71.3 kcyc   0.2% x1 119 us/call
  fx_prepare_frame       877.81 us 526.7 kcyc   1.8% x1 878 us/call
  fx_advance             922.43 us 553.5 kcyc   1.9% x1 922 us/call
    palette_cycle_step   870.01 us 522.0 kcyc   1.8% x1 870 us/call
  fx_timeline_step       111.24 us  66.7 kcyc   0.2% x1 111 us/call
  canvas_clear            86.60 us  52.0 kcyc   0.2% x1 86.6 us/call
  canvas_buffer_wait       1.99 us   1.2 kcyc   0.0% x1 1.99 us/call
```

Wall min/mean/max was 45.359/48.226/51.517 ms. Render averaged 48.225 ms per frame, including 45.778 ms in `fx_shader_draw`.

This regime leaves 14.275 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 385–416, measured over the separate 2,001,686 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1153.1/frame  0.500/1.748/17.111 us  3.22%
  isr_pack           144.1/frame  6.878/7.586/10.170 us  1.74%
  isr_dma_submit     144.1/frame  0.640/0.913/0.985 us  0.21%
```

Inclusive wake time leaves 60.487 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 46.370 ms/frame (95.1% of draw time).
2. `fx_prepare_column` took 2.499 ms/frame (5.1% of draw time).
3. `fx_advance` took 0.925 ms/frame (1.9% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Millefiori` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
