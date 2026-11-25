# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/aliencore_ship.log`. Supersedes `profile_aliencore_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | AlienCore at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh AlienCore profile 66 32` from the captured revision |

Captured 2026-11-20 15:53, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `e82fb9b061a3cb57e135d6afcc20af86410d52da399e2c3784fcf209f388d83e`.

Profile image size `FLASH: code:100648, data:152604, headers:8892 free for files:1769472 / RAM1: variables:314144, code:27352, padding:5416 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `76c4914acb13e4a29d90ef0ba11d73c76dfef6e205055f911e2073c6a7c3d6a8`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 897–928. The root counter measured 751814.233 µs, compared with 751,814 µs of draw time (0.31 ppm).

## Frame cadence

**Peak confirmed live render 27.527 ms** at frame 381 (AlienCore). **0/941 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 41.722 ms. Live render averaged 21.842 ms.

Among complete live readouts, the largest window mean was 23.492 ms/frame at frames 897–928. The dominant scope averaged 18.524 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–942 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 18 frame(s), 943–960, with uncertain membership. Their largest render was 26.730 ms and 0 exceeded the budget. The live denominator lies between 941 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 27.527, spilled 🟢 0/941 (0.00%).

README setup frame 1: 41.722 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    21.747 ms 13.05 Mcyc 100.0%
  pov_preserve_half      147.17 us  88.3 kcyc   0.7% x0.97 152 us/call
  fx_shader_draw         18.927 ms 11.36 Mcyc  87.0% x1 18927 us/call
    fx_prepare_column    615.31 us 369.2 kcyc   2.8% x144 4.27 us/call
  fx_palette_base_bake   598.88 us 359.3 kcyc   2.8% x1 599 us/call
    fx_palette_rgb_bake  484.24 us 290.5 kcyc   2.2% x1 484 us/call
  fx_prepare_frame       633.67 us 380.2 kcyc   2.9% x1 634 us/call
  fx_advance              1.029 ms 617.6 kcyc   4.7% x1 1029 us/call
    palette_cycle_step   991.26 us 594.8 kcyc   4.6% x1 991 us/call
  fx_timeline_step        80.43 us  48.3 kcyc   0.4% x1 80.4 us/call
  canvas_clear            99.75 us  59.9 kcyc   0.5% x1 99.8 us/call
  canvas_buffer_wait       2.28 us   1.4 kcyc   0.0% x1 2.28 us/call
```

Wall min/mean/max was 19.563/21.747/41.726 ms. Render averaged 21.745 ms per frame, including 18.927 ms in `fx_shader_draw`.

This regime leaves 40.755 ms of mean render headroom per display window.

### Window containing the live peak (frames 353–384)

```
frame                    22.015 ms 13.21 Mcyc 100.0%
  pov_preserve_half      145.45 us  87.3 kcyc   0.7% x1 145 us/call
  fx_shader_draw         18.461 ms 11.08 Mcyc  83.9% x1 18461 us/call
    fx_prepare_column    611.81 us 367.1 kcyc   2.8% x144 4.25 us/call
  fx_palette_base_bake    1.296 ms 777.4 kcyc   5.9% x1 1296 us/call
    fx_palette_rgb_bake   1.178 ms 706.5 kcyc   5.3% x1 1178 us/call
  fx_prepare_frame       711.62 us 427.0 kcyc   3.2% x1 712 us/call
  fx_advance             998.77 us 599.3 kcyc   4.5% x1 999 us/call
    palette_cycle_step   961.48 us 576.9 kcyc   4.4% x1 961 us/call
  fx_timeline_step        81.51 us  48.9 kcyc   0.4% x1 81.5 us/call
  canvas_clear            94.49 us  56.7 kcyc   0.4% x1 94.5 us/call
  canvas_buffer_wait       2.43 us   1.5 kcyc   0.0% x1 2.43 us/call
```

Wall min/mean/max was 20.211/22.015/27.528 ms. Render averaged 22.013 ms per frame, including 18.461 ms in `fx_shader_draw`.

This regime leaves 40.487 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    23.494 ms 14.10 Mcyc 100.0%
  pov_preserve_half      146.69 us  88.0 kcyc   0.6% x1 147 us/call
  fx_shader_draw         19.331 ms 11.60 Mcyc  82.3% x1 19331 us/call
    fx_prepare_column    614.42 us 368.7 kcyc   2.6% x144 4.27 us/call
  fx_palette_base_bake    1.953 ms  1.17 Mcyc   8.3% x1 1953 us/call
    fx_palette_rgb_bake   1.835 ms  1.10 Mcyc   7.8% x1 1835 us/call
  fx_prepare_frame       607.41 us 364.4 kcyc   2.6% x1 607 us/call
  fx_advance              1.056 ms 633.6 kcyc   4.5% x1 1056 us/call
    palette_cycle_step    1.019 ms 611.4 kcyc   4.3% x1 1019 us/call
  fx_timeline_step        78.95 us  47.4 kcyc   0.3% x1 78.9 us/call
  canvas_clear            94.03 us  56.4 kcyc   0.4% x1 94.0 us/call
  canvas_buffer_wait       2.52 us   1.5 kcyc   0.0% x1 2.52 us/call
```

Wall min/mean/max was 20.037/23.494/27.066 ms. Render averaged 23.492 ms per frame, including 19.331 ms in `fx_shader_draw`.

This regime leaves 39.008 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 353–384, measured over the separate 1,999,991 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.1/frame  0.578/1.748/16.060 us  3.22%
  isr_pack           144.0/frame  6.760/7.588/10.246 us  1.74%
  isr_dma_submit     144.0/frame  0.646/0.918/4.826 us  0.21%
```

Inclusive wake time leaves 60.487 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 18.461 ms/frame (83.9% of draw time).
2. `fx_palette_base_bake` took 1.296 ms/frame (5.9% of draw time).
3. `fx_palette_rgb_bake` took 1.178 ms/frame (5.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=AlienCore` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
