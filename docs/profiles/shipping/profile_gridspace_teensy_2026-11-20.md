# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/gridspace_ship.log`. Supersedes `profile_gridspace_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | GridSpace at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh GridSpace profile 66 32` from the captured revision |

Captured 2026-11-20 15:21, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `700ede0aefd7b0ce8c282434211f570367c832e42adaf28e2ee0abf1bc56780a`.

Profile image size `FLASH: code:99592, data:152260, headers:8244 free for files:1771520 / RAM1: variables:314144, code:27080, padding:5688 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `2d5ae33fdadbd92f32bd1e459407a2d75ed4bf50182f0c3063dd2e4ef2b9ac71`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 385–416. The root counter measured 871084.717 µs, compared with 871,084 µs of draw time (0.82 ppm).

## Frame cadence

**Peak confirmed live render 32.296 ms** at frame 403 (GridSpace). **0/943 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 45.993 ms. Live render averaged 26.120 ms.

Among complete live readouts, the largest window mean was 27.324 ms/frame at frames 897–928. The dominant scope averaged 20.439 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–944 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 945–960, with uncertain membership. Their largest render was 30.662 ms and 0 exceeded the budget. The live denominator lies between 943 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 32.296, spilled 🟢 0/943 (0.00%).

README setup frame 1: 45.993 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    25.901 ms 15.54 Mcyc 100.0%
  pov_preserve_half      139.10 us  83.5 kcyc   0.5% x0.97 144 us/call
  fx_shader_draw         20.766 ms 12.46 Mcyc  80.2% x1 20766 us/call
    fx_prepare_column    500.15 us 300.1 kcyc   1.9% x144 3.47 us/call
  fx_palette_base_bake   603.68 us 362.2 kcyc   2.3% x1 604 us/call
    fx_palette_rgb_bake  483.00 us 289.8 kcyc   1.9% x1 483 us/call
  fx_prepare_frame        2.963 ms  1.78 Mcyc  11.4% x1 2963 us/call
  fx_advance              1.069 ms 641.4 kcyc   4.1% x1 1069 us/call
    palette_cycle_step    1.001 ms 600.6 kcyc   3.9% x1 1001 us/call
  fx_timeline_step        80.94 us  48.6 kcyc   0.3% x1 80.9 us/call
  canvas_clear            96.29 us  57.8 kcyc   0.4% x1 96.3 us/call
  canvas_buffer_wait       1.80 us   1.1 kcyc   0.0% x1 1.80 us/call
```

Wall min/mean/max was 23.535/25.900/45.996 ms. Render averaged 25.899 ms per frame, including 20.766 ms in `fx_shader_draw`.

This regime leaves 36.601 ms of mean render headroom per display window.

### Window containing the live peak (frames 385–416)

```
frame                    27.221 ms 16.33 Mcyc 100.0%
  pov_preserve_half      143.32 us  86.0 kcyc   0.5% x1 143 us/call
  fx_shader_draw         21.070 ms 12.64 Mcyc  77.4% x1 21070 us/call
    fx_prepare_column    499.65 us 299.8 kcyc   1.8% x144 3.47 us/call
  fx_palette_base_bake    1.496 ms 897.3 kcyc   5.5% x1 1496 us/call
    fx_palette_rgb_bake   1.366 ms 819.7 kcyc   5.0% x1 1366 us/call
  fx_prepare_frame        3.094 ms  1.86 Mcyc  11.4% x1 3094 us/call
  fx_advance              1.074 ms 644.6 kcyc   3.9% x1 1074 us/call
    palette_cycle_step    1.003 ms 601.6 kcyc   3.7% x1 1003 us/call
  fx_timeline_step        80.79 us  48.5 kcyc   0.3% x1 80.8 us/call
  canvas_clear            92.68 us  55.6 kcyc   0.3% x1 92.7 us/call
  canvas_buffer_wait       1.72 us   1.0 kcyc   0.0% x1 1.72 us/call
```

Wall min/mean/max was 24.818/27.221/32.297 ms. Render averaged 27.220 ms per frame, including 21.070 ms in `fx_shader_draw`.

This regime leaves 35.280 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    27.326 ms 16.40 Mcyc 100.0%
  pov_preserve_half      141.69 us  85.0 kcyc   0.5% x1 142 us/call
  fx_shader_draw         20.873 ms 12.52 Mcyc  76.4% x1 20873 us/call
    fx_prepare_column    497.69 us 298.6 kcyc   1.8% x144 3.46 us/call
  fx_palette_base_bake    1.945 ms  1.17 Mcyc   7.1% x1 1945 us/call
    fx_palette_rgb_bake   1.818 ms  1.09 Mcyc   6.7% x1 1818 us/call
  fx_prepare_frame        3.014 ms  1.81 Mcyc  11.0% x1 3014 us/call
  fx_advance             999.18 us 599.5 kcyc   3.7% x1 999 us/call
    palette_cycle_step   934.64 us 560.8 kcyc   3.4% x1 935 us/call
  fx_timeline_step        85.44 us  51.3 kcyc   0.3% x1 85.4 us/call
  canvas_clear            92.84 us  55.7 kcyc   0.3% x1 92.8 us/call
  canvas_buffer_wait       1.74 us   1.0 kcyc   0.0% x1 1.74 us/call
```

Wall min/mean/max was 24.056/27.325/31.025 ms. Render averaged 27.324 ms per frame, including 20.873 ms in `fx_shader_draw`.

This regime leaves 35.176 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 385–416, measured over the separate 2,000,397 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.3/frame  0.500/1.703/24.045 us  3.14%
  isr_pack           144.0/frame  4.576/7.241/9.968 us  1.66%
  isr_dma_submit     144.0/frame  0.555/0.916/6.928 us  0.21%
```

Inclusive wake time leaves 60.538 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 21.070 ms/frame (77.4% of draw time).
2. `fx_prepare_frame` took 3.094 ms/frame (11.4% of draw time).
3. `fx_palette_base_bake` took 1.496 ms/frame (5.5% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=GridSpace` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
