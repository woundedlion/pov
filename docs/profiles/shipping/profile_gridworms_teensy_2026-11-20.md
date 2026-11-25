# GridWorms on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/gridworms_ship.log`. Supersedes `profile_gridworms_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | GridWorms at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh GridWorms profile 66 32` from the captured revision |

Captured 2026-11-20 15:59, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `3e4df00e69b1f37a101dcb2ffc9e83c264fd6874770558462b2378ddb0f59e0c`.

Profile image size `FLASH: code:98696, data:152428, headers:8972 free for files:1771520 / RAM1: variables:314144, code:26808, padding:5960 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `30b6b1b8b6c2cc3e5caee46dce9ac951c7929ae3118316d6372d431d3b3a3063`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 897–928. The root counter measured 1032665.692 µs, compared with 1,032,666 µs of draw time (0.30 ppm).

## Frame cadence

**Peak confirmed live render 38.612 ms** at frame 307 (GridWorms). **0/945 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 52.447 ms. Live render averaged 29.096 ms.

Among complete live readouts, the largest window mean was 33.107 ms/frame at frames 289–320. The dominant scope averaged 25.764 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–946 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 947–960, with uncertain membership. Their largest render was 32.668 ms and 0 exceeded the budget. The live denominator lies between 945 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 38.612, spilled 🟢 0/945 (0.00%).

README setup frame 1: 52.447 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    29.224 ms 17.53 Mcyc 100.0%
  pov_preserve_half      136.42 us  81.8 kcyc   0.5% x0.97 141 us/call
  fx_shader_draw         26.148 ms 15.69 Mcyc  89.5% x1 26148 us/call
    fx_prepare_column    666.47 us 399.9 kcyc   2.3% x144 4.63 us/call
  fx_palette_base_bake   830.73 us 498.4 kcyc   2.8% x1 831 us/call
    fx_palette_rgb_bake  702.17 us 421.3 kcyc   2.4% x1 702 us/call
  fx_prepare_frame       727.86 us 436.7 kcyc   2.5% x1 728 us/call
  fx_advance             952.53 us 571.5 kcyc   3.3% x1 953 us/call
    palette_cycle_step   898.22 us 538.9 kcyc   3.1% x1 898 us/call
  fx_timeline_step       112.39 us  67.4 kcyc   0.4% x1 112 us/call
  canvas_clear            99.83 us  59.9 kcyc   0.3% x1 99.8 us/call
  canvas_buffer_wait       1.62 us    972 cyc   0.0% x1 1.62 us/call
```

Wall min/mean/max was 25.357/29.223/52.448 ms. Render averaged 29.222 ms per frame, including 26.148 ms in `fx_shader_draw`.

This regime leaves 33.278 ms of mean render headroom per display window.

### Window containing the live peak (frames 289–320)

```
frame                    33.109 ms 19.87 Mcyc 100.0%
  pov_preserve_half      141.43 us  84.9 kcyc   0.4% x1 141 us/call
  fx_shader_draw         28.534 ms 17.12 Mcyc  86.2% x1 28534 us/call
    fx_prepare_column    775.51 us 465.3 kcyc   2.3% x144 5.39 us/call
  fx_palette_base_bake    2.318 ms  1.39 Mcyc   7.0% x1 2318 us/call
    fx_palette_rgb_bake   2.191 ms  1.31 Mcyc   6.6% x1 2191 us/call
  fx_prepare_frame       793.46 us 476.1 kcyc   2.4% x1 793 us/call
  fx_advance             917.06 us 550.2 kcyc   2.8% x1 917 us/call
    palette_cycle_step   864.08 us 518.4 kcyc   2.6% x1 864 us/call
  fx_timeline_step       101.02 us  60.6 kcyc   0.3% x1 101 us/call
  canvas_clear            96.06 us  57.6 kcyc   0.3% x1 96.1 us/call
  canvas_buffer_wait       1.68 us   1.0 kcyc   0.0% x1 1.68 us/call
```

Wall min/mean/max was 28.997/33.108/38.614 ms. Render averaged 33.107 ms per frame, including 28.534 ms in `fx_shader_draw`.

This regime leaves 29.393 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    32.271 ms 19.36 Mcyc 100.0%
  pov_preserve_half      142.27 us  85.4 kcyc   0.4% x1 142 us/call
  fx_shader_draw         28.868 ms 17.32 Mcyc  89.5% x1 28868 us/call
    fx_prepare_column    787.57 us 472.5 kcyc   2.4% x144 5.47 us/call
  fx_palette_base_bake    1.251 ms 750.6 kcyc   3.9% x1 1251 us/call
    fx_palette_rgb_bake   1.119 ms 671.2 kcyc   3.5% x1 1119 us/call
  fx_prepare_frame       674.96 us 405.0 kcyc   2.1% x1 675 us/call
  fx_advance             916.67 us 550.0 kcyc   2.8% x1 917 us/call
    palette_cycle_step   862.72 us 517.6 kcyc   2.7% x1 863 us/call
  fx_timeline_step       100.61 us  60.4 kcyc   0.3% x1 101 us/call
  canvas_clear            96.83 us  58.1 kcyc   0.3% x1 96.8 us/call
  canvas_buffer_wait       1.58 us    947 cyc   0.0% x1 1.58 us/call
```

Wall min/mean/max was 28.868/32.270/33.425 ms. Render averaged 32.269 ms per frame, including 28.868 ms in `fx_shader_draw`.

This regime leaves 30.231 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 289–320, measured over the separate 1,999,462 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.8/frame  0.521/1.883/21.805 us  3.47%
  isr_pack           144.0/frame  5.268/7.191/10.021 us  1.65%
  isr_dma_submit     144.0/frame  0.605/0.913/1.070 us  0.21%
```

Inclusive wake time leaves 60.331 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 28.534 ms/frame (86.2% of draw time).
2. `fx_palette_base_bake` took 2.318 ms/frame (7.0% of draw time).
3. `fx_palette_rgb_bake` took 2.191 ms/frame (6.6% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=GridWorms` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
