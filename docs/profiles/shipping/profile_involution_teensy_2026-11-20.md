# Involution on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/involution_ship.log`. Supersedes `profile_involution_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Involution at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Involution profile 66 32` from the captured revision |

Captured 2026-11-20 16:00, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `056a79befb5c045deadff2378843768f4ceaa34525c22e62de7ac8e341233c4a`.

Profile image size `FLASH: code:113264, data:159804, headers:8532 free for files:1750016 / RAM1: variables:314144, code:27336, padding:5432 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `4fafcd48e62123b7c7522c8165b874cf2b7428949eb4741d0a0ec62fff2c60ab`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 33–64. The root counter measured 1433835.488 µs, compared with 1,433,835 µs of draw time (0.34 ppm).

## Frame cadence

**Peak confirmed live render 45.232 ms** at frame 46 (Involution). **0/949 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 81.961 ms. Live render averaged 44.537 ms.

Among complete live readouts, the largest window mean was 44.805 ms/frame at frames 33–64. The dominant scope averaged 43.927 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–950 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 10 frame(s), 951–960, with uncertain membership. Their largest render was 44.455 ms and 0 exceeded the budget. The live denominator lies between 949 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 45.232, spilled 🟢 0/949 (0.00%).

README setup frame 1: 81.961 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    45.960 ms 27.58 Mcyc 100.0%
  pov_preserve_half      132.71 us  79.6 kcyc   0.3% x0.97 137 us/call
  fx_shader_draw         45.344 ms 27.21 Mcyc  98.7% x1 45344 us/call
    fx_prepare_column     1.454 ms 872.2 kcyc   3.2% x144 10.1 us/call
  fx_palette_base_bake    69.69 us  41.8 kcyc   0.2% x1 69.7 us/call
  fx_prepare_frame        16.74 us  10.0 kcyc   0.0% x1 16.7 us/call
  fx_advance              50.69 us  30.4 kcyc   0.1% x1 50.7 us/call
  fx_timeline_step        49.37 us  29.6 kcyc   0.1% x1 49.4 us/call
  canvas_clear            90.95 us  54.6 kcyc   0.2% x1 91.0 us/call
  canvas_buffer_wait       2.76 us   1.7 kcyc   0.0% x1 2.76 us/call
```

Wall min/mean/max was 43.991/45.959/81.961 ms. Render averaged 45.957 ms per frame, including 45.344 ms in `fx_shader_draw`.

This regime leaves 16.543 ms of mean render headroom per display window.

### Window containing the live peak (frames 33–64)

```
frame                    44.807 ms 26.88 Mcyc 100.0%
  pov_preserve_half      138.53 us  83.1 kcyc   0.3% x1 139 us/call
  fx_shader_draw         44.197 ms 26.52 Mcyc  98.6% x1 44197 us/call
    fx_prepare_column     1.523 ms 913.7 kcyc   3.4% x144 10.6 us/call
  fx_palette_base_bake    68.64 us  41.2 kcyc   0.2% x1 68.6 us/call
  fx_prepare_frame        15.63 us   9.4 kcyc   0.0% x1 15.6 us/call
  fx_advance              47.16 us  28.3 kcyc   0.1% x1 47.2 us/call
  fx_timeline_step        48.66 us  29.2 kcyc   0.1% x1 48.7 us/call
  canvas_clear            86.47 us  51.9 kcyc   0.2% x1 86.5 us/call
  canvas_buffer_wait       1.91 us   1.1 kcyc   0.0% x1 1.91 us/call
```

Wall min/mean/max was 44.341/44.807/45.233 ms. Render averaged 44.805 ms per frame, including 44.197 ms in `fx_shader_draw`.

This regime leaves 17.695 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    44.534 ms 26.72 Mcyc 100.0%
  pov_preserve_half      137.19 us  82.3 kcyc   0.3% x1 137 us/call
  fx_shader_draw         43.924 ms 26.35 Mcyc  98.6% x1 43924 us/call
    fx_prepare_column     1.589 ms 953.7 kcyc   3.6% x144 11.0 us/call
  fx_palette_base_bake    68.07 us  40.8 kcyc   0.2% x1 68.1 us/call
  fx_prepare_frame        15.49 us   9.3 kcyc   0.0% x1 15.5 us/call
  fx_advance              48.66 us  29.2 kcyc   0.1% x1 48.7 us/call
  fx_timeline_step        47.38 us  28.4 kcyc   0.1% x1 47.4 us/call
  canvas_clear            86.66 us  52.0 kcyc   0.2% x1 86.7 us/call
  canvas_buffer_wait       2.73 us   1.6 kcyc   0.0% x1 2.73 us/call
```

Wall min/mean/max was 44.254/44.533/44.890 ms. Render averaged 44.531 ms per frame, including 43.924 ms in `fx_shader_draw`.

This regime leaves 17.969 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 33–64, measured over the separate 1,999,871 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.0/frame  0.500/1.943/17.931 us  3.58%
  isr_pack           144.0/frame  6.416/7.535/10.078 us  1.73%
  isr_dma_submit     144.0/frame  0.615/0.915/9.483 us  0.21%
```

Inclusive wake time leaves 60.262 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 44.197 ms/frame (98.6% of draw time).
2. `fx_prepare_column` took 1.523 ms/frame (3.4% of draw time).
3. `pov_preserve_half` took 0.139 ms/frame (0.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Involution` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
