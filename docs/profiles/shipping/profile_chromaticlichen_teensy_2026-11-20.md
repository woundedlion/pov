# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/chromaticlichen_ship.log`. Supersedes `profile_chromaticlichen_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | ChromaticLichen at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh ChromaticLichen profile 66 32` from the captured revision |

Captured 2026-11-20 16:10, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `6c4b75cd932dac00c38c91b47147d9436173f81b48547b0a038705d537fe7cf2`.

Profile image size `FLASH: code:101376, data:151644, headers:9124 free for files:1769472 / RAM1: variables:314144, code:28760, padding:4008 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `9ada6cc1b8648b8d0eb4e954e2f832b87b3535a3f6a8035404edd43a4b8360a3`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 321–352. The root counter measured 1326161.267 µs, compared with 1,326,161 µs of draw time (0.20 ppm).

## Frame cadence

**Peak confirmed live render 45.706 ms** at frame 259 (ChromaticLichen). **0/948 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 68.909 ms. Live render averaged 37.928 ms.

Among complete live readouts, the largest window mean was 41.437 ms/frame at frames 321–352. The dominant scope averaged 32.226 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–949 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 11 frame(s), 950–960, with uncertain membership. Their largest render was 40.087 ms and 0 exceeded the budget. The live denominator lies between 948 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 45.706, spilled 🟢 0/948 (0.00%).

README setup frame 1: 68.909 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    38.558 ms 23.13 Mcyc 100.0%
  pov_preserve_half      134.06 us  80.4 kcyc   0.3% x0.97 138 us/call
  fx_shader_draw         33.126 ms 19.88 Mcyc  85.9% x1 33126 us/call
    fx_prepare_column    561.05 us 336.6 kcyc   1.5% x144 3.90 us/call
  fx_palette_base_bake   814.44 us 488.7 kcyc   2.1% x1 814 us/call
    fx_palette_rgb_bake  695.96 us 417.6 kcyc   1.8% x1 696 us/call
  fx_prepare_frame        3.117 ms  1.87 Mcyc   8.1% x1 3117 us/call
  fx_advance             975.89 us 585.5 kcyc   2.5% x1 976 us/call
    palette_cycle_step   907.85 us 544.7 kcyc   2.4% x1 908 us/call
  fx_timeline_step        95.69 us  57.4 kcyc   0.2% x1 95.7 us/call
  canvas_clear            99.07 us  59.4 kcyc   0.3% x1 99.1 us/call
  canvas_buffer_wait       2.62 us   1.6 kcyc   0.0% x1 2.62 us/call
```

Wall min/mean/max was 35.770/38.558/68.911 ms. Render averaged 38.555 ms per frame, including 33.126 ms in `fx_shader_draw`.

This regime leaves 23.945 ms of mean render headroom per display window.

### Window containing the live peak (frames 257–288)

```
frame                    40.484 ms 24.29 Mcyc 100.0%
  pov_preserve_half      140.72 us  84.4 kcyc   0.3% x1 141 us/call
  fx_shader_draw         33.984 ms 20.39 Mcyc  83.9% x1 33984 us/call
    fx_prepare_column    561.64 us 337.0 kcyc   1.4% x144 3.90 us/call
  fx_palette_base_bake    1.807 ms  1.08 Mcyc   4.5% x1 1807 us/call
    fx_palette_rgb_bake   1.689 ms  1.01 Mcyc   4.2% x1 1689 us/call
  fx_prepare_frame        3.243 ms  1.95 Mcyc   8.0% x1 3243 us/call
  fx_advance             919.61 us 551.8 kcyc   2.3% x1 920 us/call
    palette_cycle_step   851.20 us 510.7 kcyc   2.1% x1 851 us/call
  fx_timeline_step        95.47 us  57.3 kcyc   0.2% x1 95.5 us/call
  canvas_clear            94.09 us  56.5 kcyc   0.2% x1 94.1 us/call
  canvas_buffer_wait       4.05 us   2.4 kcyc   0.0% x1 4.05 us/call
```

Wall min/mean/max was 38.311/40.483/45.710 ms. Render averaged 40.480 ms per frame, including 33.984 ms in `fx_shader_draw`.

This regime leaves 22.020 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    38.289 ms 22.97 Mcyc 100.0%
  pov_preserve_half      139.02 us  83.4 kcyc   0.4% x1 139 us/call
  fx_shader_draw         32.469 ms 19.48 Mcyc  84.8% x1 32469 us/call
    fx_prepare_column    571.16 us 342.7 kcyc   1.5% x144 3.97 us/call
  fx_palette_base_bake    1.236 ms 741.5 kcyc   3.2% x1 1236 us/call
    fx_palette_rgb_bake   1.120 ms 671.8 kcyc   2.9% x1 1120 us/call
  fx_prepare_frame        3.125 ms  1.87 Mcyc   8.2% x1 3125 us/call
  fx_advance             929.27 us 557.6 kcyc   2.4% x1 929 us/call
    palette_cycle_step   863.22 us 517.9 kcyc   2.3% x1 863 us/call
  fx_timeline_step        97.06 us  58.2 kcyc   0.3% x1 97.1 us/call
  canvas_clear            94.75 us  56.9 kcyc   0.2% x1 94.8 us/call
  canvas_buffer_wait       5.53 us   3.3 kcyc   0.0% x1 5.53 us/call
```

Wall min/mean/max was 36.647/38.289/39.844 ms. Render averaged 38.284 ms per frame, including 32.469 ms in `fx_shader_draw`.

This regime leaves 24.216 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 257–288, measured over the separate 1,999,751 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.9/frame  0.500/1.711/21.858 us  3.15%
  isr_pack           144.0/frame  6.401/7.208/9.965 us  1.66%
  isr_dma_submit     144.0/frame  0.593/0.915/6.961 us  0.21%
```

Inclusive wake time leaves 60.531 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 33.984 ms/frame (83.9% of draw time).
2. `fx_prepare_frame` took 3.243 ms/frame (8.0% of draw time).
3. `fx_palette_base_bake` took 1.807 ms/frame (4.5% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=ChromaticLichen` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
