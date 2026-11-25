# JewelMelt on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/jewelmelt_ship.log`. Supersedes `profile_jewelmelt_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | JewelMelt at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh JewelMelt profile 66 32` from the captured revision |

Captured 2026-11-20 16:01, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `eafc9b0f7d920b11515ca32dd918b0d3c3f5feda0e6c8eae0ad04c84aea2abdf`.

Profile image size `FLASH: code:101464, data:151720, headers:8960 free for files:1769472 / RAM1: variables:314144, code:28760, padding:4008 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `c469b5b588d53d90978c5f9d9abc78684dc020b6edbccc27345c88ab561b8028`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 481–512. The root counter measured 1486963.025 µs, compared with 1,486,965 µs of draw time (1.33 ppm).

## Frame cadence

**Peak confirmed live render 53.778 ms** at frame 319 (JewelMelt). **0/949 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 83.605 ms. Live render averaged 45.185 ms.

Among complete live readouts, the largest window mean was 48.417 ms/frame at frames 321–352. The dominant scope averaged 39.264 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–950 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 10 frame(s), 951–960, with uncertain membership. Their largest render was 47.352 ms and 0 exceeded the budget. The live denominator lies between 949 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 53.778, spilled 🟢 0/949 (0.00%).

README setup frame 1: 83.605 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    45.678 ms 27.41 Mcyc 100.0%
  pov_preserve_half      136.72 us  82.0 kcyc   0.3% x0.97 141 us/call
  fx_shader_draw         39.886 ms 23.93 Mcyc  87.3% x1 39886 us/call
    fx_prepare_column    506.56 us 303.9 kcyc   1.1% x144 3.52 us/call
  fx_palette_base_bake   811.23 us 486.7 kcyc   1.8% x1 811 us/call
    fx_palette_rgb_bake  703.05 us 421.8 kcyc   1.5% x1 703 us/call
  fx_prepare_frame        3.165 ms  1.90 Mcyc   6.9% x1 3165 us/call
  fx_advance              1.320 ms 792.1 kcyc   2.9% x1 1320 us/call
    palette_cycle_step    1.255 ms 752.9 kcyc   2.7% x1 1255 us/call
  fx_timeline_step        83.57 us  50.1 kcyc   0.2% x1 83.6 us/call
  canvas_clear            99.46 us  59.7 kcyc   0.2% x1 99.5 us/call
  canvas_buffer_wait       1.53 us    919 cyc   0.0% x1 1.53 us/call
```

Wall min/mean/max was 42.135/45.678/83.607 ms. Render averaged 45.676 ms per frame, including 39.886 ms in `fx_shader_draw`.

This regime leaves 16.824 ms of mean render headroom per display window.

### Window containing the live peak (frames 289–320)

```
frame                    47.290 ms 28.37 Mcyc 100.0%
  pov_preserve_half      139.63 us  83.8 kcyc   0.3% x1 140 us/call
  fx_shader_draw         40.095 ms 24.06 Mcyc  84.8% x1 40095 us/call
    fx_prepare_column    507.40 us 304.4 kcyc   1.1% x144 3.52 us/call
  fx_palette_base_bake    2.295 ms  1.38 Mcyc   4.9% x1 2295 us/call
    fx_palette_rgb_bake   2.183 ms  1.31 Mcyc   4.6% x1 2183 us/call
  fx_prepare_frame        3.223 ms  1.93 Mcyc   6.8% x1 3223 us/call
  fx_advance              1.195 ms 716.8 kcyc   2.5% x1 1195 us/call
    palette_cycle_step    1.129 ms 677.2 kcyc   2.4% x1 1129 us/call
  fx_timeline_step        72.68 us  43.6 kcyc   0.2% x1 72.7 us/call
  canvas_clear            91.66 us  55.0 kcyc   0.2% x1 91.7 us/call
  canvas_buffer_wait       1.26 us    754 cyc   0.0% x1 1.26 us/call
```

Wall min/mean/max was 42.535/47.290/53.779 ms. Render averaged 47.289 ms per frame, including 40.095 ms in `fx_shader_draw`.

This regime leaves 15.211 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    46.375 ms 27.82 Mcyc 100.0%
  pov_preserve_half      138.90 us  83.3 kcyc   0.3% x1 139 us/call
  fx_shader_draw         40.352 ms 24.21 Mcyc  87.0% x1 40352 us/call
    fx_prepare_column    540.64 us 324.4 kcyc   1.2% x144 3.75 us/call
  fx_palette_base_bake    1.230 ms 738.2 kcyc   2.7% x1 1230 us/call
    fx_palette_rgb_bake   1.118 ms 671.1 kcyc   2.4% x1 1118 us/call
  fx_prepare_frame        3.129 ms  1.88 Mcyc   6.7% x1 3129 us/call
  fx_advance              1.189 ms 713.2 kcyc   2.6% x1 1189 us/call
    palette_cycle_step    1.123 ms 673.8 kcyc   2.4% x1 1123 us/call
  fx_timeline_step        73.34 us  44.0 kcyc   0.2% x1 73.3 us/call
  canvas_clear            91.21 us  54.7 kcyc   0.2% x1 91.2 us/call
  canvas_buffer_wait       1.66 us    997 cyc   0.0% x1 1.66 us/call
```

Wall min/mean/max was 42.615/46.375/50.245 ms. Render averaged 46.373 ms per frame, including 40.352 ms in `fx_shader_draw`.

This regime leaves 16.127 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 289–320, measured over the separate 2,001,419 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.9/frame  0.578/1.758/28.653 us  3.24%
  isr_pack           144.1/frame  6.845/7.555/10.453 us  1.74%
  isr_dma_submit     144.1/frame  0.628/0.916/1.205 us  0.21%
```

Inclusive wake time leaves 60.475 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 40.095 ms/frame (84.8% of draw time).
2. `fx_prepare_frame` took 3.223 ms/frame (6.8% of draw time).
3. `fx_palette_base_bake` took 2.295 ms/frame (4.9% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=JewelMelt` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
