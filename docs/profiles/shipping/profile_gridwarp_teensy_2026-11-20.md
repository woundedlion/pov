# GridWarp on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/gridwarp_ship.log`. Supersedes `profile_gridwarp_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | GridWarp at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh GridWarp profile 66 32` from the captured revision |

Captured 2026-11-20 16:14, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `aca3af2a84df08251bf5b0e6d40953f41f94d2105e637368385f1ac7ecb6c149`.

Profile image size `FLASH: code:110600, data:154908, headers:8924 free for files:1757184 / RAM1: variables:314144, code:27304, padding:5464 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `f89f721b9d7beba8998f82efbb10ebe176cd71b7f7debc582e167febced790ef`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 225–256. The root counter measured 1459448.658 µs, compared with 1,459,454 µs of draw time (3.66 ppm).

## Frame cadence

**Peak confirmed live render 50.732 ms** at frame 247 (GridWarp). **0/949 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 78.390 ms. Live render averaged 42.704 ms.

Among complete live readouts, the largest window mean was 45.606 ms/frame at frames 225–256. The dominant scope averaged 36.739 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–950 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 10 frame(s), 951–960, with uncertain membership. Their largest render was 46.079 ms and 0 exceeded the budget. The live denominator lies between 949 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 50.732, spilled 🟢 0/949 (0.00%).

README setup frame 1: 78.390 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    44.046 ms 26.43 Mcyc 100.0%
  pov_preserve_half      135.34 us  81.2 kcyc   0.3% x0.97 140 us/call
  fx_shader_draw         38.232 ms 22.94 Mcyc  86.8% x1 38232 us/call
    fx_prepare_column     1.205 ms 723.2 kcyc   2.7% x144 8.37 us/call
  fx_palette_base_bake   831.06 us 498.6 kcyc   1.9% x1 831 us/call
    fx_palette_rgb_bake  708.76 us 425.3 kcyc   1.6% x1 709 us/call
  fx_prepare_frame        3.433 ms  2.06 Mcyc   7.8% x1 3433 us/call
  fx_advance              1.078 ms 646.8 kcyc   2.4% x1 1078 us/call
    palette_cycle_step    1.038 ms 622.7 kcyc   2.4% x1 1038 us/call
  fx_timeline_step        83.05 us  49.8 kcyc   0.2% x1 83.1 us/call
  canvas_clear            98.31 us  59.0 kcyc   0.2% x1 98.3 us/call
  canvas_buffer_wait       2.14 us   1.3 kcyc   0.0% x1 2.14 us/call
```

Wall min/mean/max was 40.382/44.045/78.395 ms. Render averaged 44.043 ms per frame, including 38.232 ms in `fx_shader_draw`.

This regime leaves 18.457 ms of mean render headroom per display window.

### Window containing the live peak (frames 225–256)

```
frame                    45.608 ms 27.36 Mcyc 100.0%
  pov_preserve_half      139.46 us  83.7 kcyc   0.3% x1 139 us/call
  fx_shader_draw         38.595 ms 23.16 Mcyc  84.6% x1 38595 us/call
    fx_prepare_column     1.208 ms 724.9 kcyc   2.6% x144 8.39 us/call
  fx_palette_base_bake    2.087 ms  1.25 Mcyc   4.6% x1 2087 us/call
    fx_palette_rgb_bake   1.962 ms  1.18 Mcyc   4.3% x1 1962 us/call
  fx_prepare_frame        3.571 ms  2.14 Mcyc   7.8% x1 3571 us/call
  fx_advance             877.26 us 526.4 kcyc   1.9% x1 877 us/call
    palette_cycle_step   836.35 us 501.8 kcyc   1.8% x1 836 us/call
  fx_timeline_step        85.43 us  51.3 kcyc   0.2% x1 85.4 us/call
  canvas_clear            93.51 us  56.1 kcyc   0.2% x1 93.5 us/call
  canvas_buffer_wait       1.96 us   1.2 kcyc   0.0% x1 1.96 us/call
```

Wall min/mean/max was 42.786/45.607/50.734 ms. Render averaged 45.606 ms per frame, including 38.595 ms in `fx_shader_draw`.

This regime leaves 16.894 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    43.461 ms 26.08 Mcyc 100.0%
  pov_preserve_half      139.77 us  83.9 kcyc   0.3% x1 140 us/call
  fx_shader_draw         37.419 ms 22.45 Mcyc  86.1% x1 37419 us/call
    fx_prepare_column     1.208 ms 724.9 kcyc   2.8% x144 8.39 us/call
  fx_palette_base_bake    1.272 ms 763.4 kcyc   2.9% x1 1272 us/call
    fx_palette_rgb_bake   1.147 ms 688.0 kcyc   2.6% x1 1147 us/call
  fx_prepare_frame        3.375 ms  2.03 Mcyc   7.8% x1 3375 us/call
  fx_advance             923.88 us 554.3 kcyc   2.1% x1 924 us/call
    palette_cycle_step   884.12 us 530.5 kcyc   2.0% x1 884 us/call
  fx_timeline_step        83.05 us  49.8 kcyc   0.2% x1 83.1 us/call
  canvas_clear            94.22 us  56.5 kcyc   0.2% x1 94.2 us/call
  canvas_buffer_wait       2.16 us   1.3 kcyc   0.0% x1 2.16 us/call
```

Wall min/mean/max was 41.283/43.461/45.687 ms. Render averaged 43.459 ms per frame, including 37.419 ms in `fx_shader_draw`.

This regime leaves 19.041 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 225–256, measured over the separate 1,999,956 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.1/frame  0.580/1.756/23.446 us  3.23%
  isr_pack           144.0/frame  6.805/7.470/10.118 us  1.72%
  isr_dma_submit     144.0/frame  0.603/0.916/4.195 us  0.21%
```

Inclusive wake time leaves 60.481 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 38.595 ms/frame (84.6% of draw time).
2. `fx_prepare_frame` took 3.571 ms/frame (7.8% of draw time).
3. `fx_palette_base_bake` took 2.087 ms/frame (4.6% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=GridWarp` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
