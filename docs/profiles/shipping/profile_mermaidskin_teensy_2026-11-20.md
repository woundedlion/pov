# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/mermaidskin_ship.log`. Supersedes `profile_mermaidskin_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | MermaidSkin at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh MermaidSkin profile 66 32` from the captured revision |

Captured 2026-11-20 16:19, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `0084fd2e3258c65bb082e4c78d41ac8c9c139e0de0b289395fb41372f7099f2a`.

Profile image size `FLASH: code:101568, data:151640, headers:8936 free for files:1769472 / RAM1: variables:314144, code:28760, padding:4008 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `f7c3e1ab4b8abcaebac979c9bb17edea925da99533ef62ad7155fb610137c045`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 225–256. The root counter measured 1311843.040 µs, compared with 1,311,842 µs of draw time (0.79 ppm).

## Frame cadence

**Peak confirmed live render 45.788 ms** at frame 249 (MermaidSkin). **0/948 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 69.784 ms. Live render averaged 38.020 ms.

Among complete live readouts, the largest window mean was 40.993 ms/frame at frames 225–256. The dominant scope averaged 32.270 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–949 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 11 frame(s), 950–960, with uncertain membership. Their largest render was 40.797 ms and 0 exceeded the budget. The live denominator lies between 948 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 45.788, spilled 🟢 0/948 (0.00%).

README setup frame 1: 69.784 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    38.949 ms 23.37 Mcyc 100.0%
  pov_preserve_half      136.66 us  82.0 kcyc   0.4% x0.97 141 us/call
  fx_shader_draw         33.479 ms 20.09 Mcyc  86.0% x1 33479 us/call
    fx_prepare_column    562.92 us 337.8 kcyc   1.4% x144 3.91 us/call
  fx_palette_base_bake   818.70 us 491.2 kcyc   2.1% x1 819 us/call
    fx_palette_rgb_bake  701.24 us 420.7 kcyc   1.8% x1 701 us/call
  fx_prepare_frame        3.134 ms  1.88 Mcyc   8.0% x1 3134 us/call
  fx_advance             978.89 us 587.3 kcyc   2.5% x1 979 us/call
    palette_cycle_step   916.13 us 549.7 kcyc   2.4% x1 916 us/call
  fx_timeline_step        86.79 us  52.1 kcyc   0.2% x1 86.8 us/call
  canvas_clear            99.18 us  59.5 kcyc   0.3% x1 99.2 us/call
  canvas_buffer_wait       1.52 us    911 cyc   0.0% x1 1.52 us/call
```

Wall min/mean/max was 36.222/38.948/69.787 ms. Render averaged 38.947 ms per frame, including 33.479 ms in `fx_shader_draw`.

This regime leaves 23.553 ms of mean render headroom per display window.

### Window containing the live peak (frames 225–256)

```
frame                    40.995 ms 24.60 Mcyc 100.0%
  pov_preserve_half      141.75 us  85.0 kcyc   0.3% x1 142 us/call
  fx_shader_draw         34.243 ms 20.55 Mcyc  83.5% x1 34243 us/call
    fx_prepare_column    558.19 us 334.9 kcyc   1.4% x144 3.88 us/call
  fx_palette_base_bake    2.056 ms  1.23 Mcyc   5.0% x1 2056 us/call
    fx_palette_rgb_bake   1.939 ms  1.16 Mcyc   4.7% x1 1939 us/call
  fx_prepare_frame        3.255 ms  1.95 Mcyc   7.9% x1 3255 us/call
  fx_advance             911.22 us 546.7 kcyc   2.2% x1 911 us/call
    palette_cycle_step   845.16 us 507.1 kcyc   2.1% x1 845 us/call
  fx_timeline_step        79.98 us  48.0 kcyc   0.2% x1 80.0 us/call
  canvas_clear            95.25 us  57.1 kcyc   0.2% x1 95.2 us/call
  canvas_buffer_wait       2.25 us   1.4 kcyc   0.0% x1 2.25 us/call
```

Wall min/mean/max was 38.585/40.995/45.789 ms. Render averaged 40.993 ms per frame, including 34.243 ms in `fx_shader_draw`.

This regime leaves 21.507 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    38.425 ms 23.05 Mcyc 100.0%
  pov_preserve_half      144.06 us  86.4 kcyc   0.4% x1 144 us/call
  fx_shader_draw         32.582 ms 19.55 Mcyc  84.8% x1 32582 us/call
    fx_prepare_column    564.14 us 338.5 kcyc   1.5% x144 3.92 us/call
  fx_palette_base_bake    1.235 ms 740.8 kcyc   3.2% x1 1235 us/call
    fx_palette_rgb_bake   1.119 ms 671.5 kcyc   2.9% x1 1119 us/call
  fx_prepare_frame        3.143 ms  1.89 Mcyc   8.2% x1 3143 us/call
  fx_advance             935.06 us 561.0 kcyc   2.4% x1 935 us/call
    palette_cycle_step   868.53 us 521.1 kcyc   2.3% x1 869 us/call
  fx_timeline_step        78.65 us  47.2 kcyc   0.2% x1 78.6 us/call
  canvas_clear            95.21 us  57.1 kcyc   0.2% x1 95.2 us/call
  canvas_buffer_wait       1.83 us   1.1 kcyc   0.0% x1 1.83 us/call
```

Wall min/mean/max was 36.473/38.424/40.080 ms. Render averaged 38.423 ms per frame, including 32.582 ms in `fx_shader_draw`.

This regime leaves 24.077 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 225–256, measured over the separate 1,999,901 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.0/frame  0.500/1.768/17.691 us  3.26%
  isr_pack           144.0/frame  6.930/7.708/10.580 us  1.77%
  isr_dma_submit     144.0/frame  0.633/0.915/4.788 us  0.21%
```

Inclusive wake time leaves 60.463 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 34.243 ms/frame (83.5% of draw time).
2. `fx_prepare_frame` took 3.255 ms/frame (7.9% of draw time).
3. `fx_palette_base_bake` took 2.056 ms/frame (5.0% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=MermaidSkin` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
