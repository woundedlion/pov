# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/ashcloud_ship.log`. Supersedes `profile_ashcloud_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | AshCloud at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh AshCloud profile 66 32` from the captured revision |

Captured 2026-11-20 14:53, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `3902934e3b8423ad7beb9a1453e66d81b8c98092e78f93f5bf25125cc305170d`.

Profile image size `FLASH: code:124568, data:163164, headers:8204 free for files:1735680 / RAM1: variables:314144, code:31704, padding:1064 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `4e93f907b8dcd9c6fed9582fcbbae8e24e4f073f0c78a96d80e06cf690c8dcb0`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 865–896. The root counter measured 1666859.152 µs, compared with 1,666,856 µs of draw time (1.89 ppm).

## Frame cadence

**Peak confirmed live render 55.513 ms** at frame 906 (AshCloud). **0/950 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 83.423 ms. Live render averaged 48.977 ms.

Among complete live readouts, the largest window mean was 52.087 ms/frame at frames 865–896. The dominant scope averaged 45.441 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–951 are confirmed live. The raw file contains 1,032 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 8 frame(s), 952–959, with uncertain membership. Their largest render was 51.295 ms and 0 exceeded the budget. The live denominator lies between 950 and 958. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 55.513, spilled 🟢 0/950 (0.00%).

README setup frame 1: 83.423 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    44.046 ms 26.43 Mcyc 100.0%
  pov_preserve_half      135.98 us  81.6 kcyc   0.3% x0.97 140 us/call
  fx_shader_draw         40.662 ms 24.40 Mcyc  92.3% x1 40662 us/call
    fx_prepare_column     1.942 ms  1.17 Mcyc   4.4% x144 13.5 us/call
  fx_palette_base_bake   589.88 us 353.9 kcyc   1.3% x1 590 us/call
    fx_palette_rgb_bake  447.38 us 268.4 kcyc   1.0% x1 447 us/call
  fx_prepare_frame        1.139 ms 683.2 kcyc   2.6% x1 1139 us/call
  fx_advance              1.140 ms 683.9 kcyc   2.6% x1 1140 us/call
    palette_cycle_step    1.068 ms 640.7 kcyc   2.4% x1 1068 us/call
  fx_timeline_step        99.19 us  59.5 kcyc   0.2% x1 99.2 us/call
  canvas_clear            90.42 us  54.2 kcyc   0.2% x1 90.4 us/call
  canvas_buffer_wait       1.96 us   1.2 kcyc   0.0% x1 1.96 us/call
```

Wall min/mean/max was 37.942/44.045/83.424 ms. Render averaged 44.044 ms per frame, including 40.662 ms in `fx_shader_draw`.

This regime leaves 18.456 ms of mean render headroom per display window.

### Window containing the live peak (frames 897–928)

```
frame                    52.069 ms 31.24 Mcyc 100.0%
  pov_preserve_half      136.40 us  81.8 kcyc   0.3% x1 136 us/call
  fx_shader_draw         48.164 ms 28.90 Mcyc  92.5% x1 48164 us/call
    fx_prepare_column     2.112 ms  1.27 Mcyc   4.1% x144 14.7 us/call
  fx_palette_base_bake   964.38 us 578.6 kcyc   1.9% x1 964 us/call
    fx_palette_rgb_bake  802.12 us 481.3 kcyc   1.5% x1 802 us/call
  fx_prepare_frame        1.325 ms 794.9 kcyc   2.5% x1 1325 us/call
  fx_advance              1.083 ms 649.9 kcyc   2.1% x1 1083 us/call
    palette_cycle_step    1.004 ms 602.7 kcyc   1.9% x1 1004 us/call
  fx_timeline_step       105.72 us  63.4 kcyc   0.2% x1 106 us/call
  canvas_clear            85.93 us  51.6 kcyc   0.2% x1 85.9 us/call
  canvas_buffer_wait       2.12 us   1.3 kcyc   0.0% x1 2.12 us/call
```

Wall min/mean/max was 49.181/52.069/55.514 ms. Render averaged 52.067 ms per frame, including 48.164 ms in `fx_shader_draw`.

This regime leaves 10.433 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 897–928, measured over the separate 2,001,304 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.8/frame  0.533/1.820/20.680 us  3.35%
  isr_pack           144.1/frame  7.031/7.728/10.625 us  1.78%
  isr_dma_submit     144.1/frame  0.628/0.925/8.410 us  0.21%
```

Inclusive wake time leaves 60.406 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 48.164 ms/frame (92.5% of draw time).
2. `fx_prepare_column` took 2.112 ms/frame (4.1% of draw time).
3. `fx_prepare_frame` took 1.325 ms/frame (2.5% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=AshCloud` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
