# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/gsreactiondiffusion_ship.log`. Supersedes `profile_gsreactiondiffusion_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | GSReactionDiffusion at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh GSReactionDiffusion profile 66 32` from the captured revision |

Captured 2026-11-20 15:56, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `4c490edc954296d7656351935a11abeb64fa0ad2ffc37c7a81610a5de4243a15`.

Profile image size `FLASH: code:80760, data:249096, headers:9088 free for files:1692672 / RAM1: variables:314240, code:18072, padding:14696 free for local variables:177280 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `385707837ae2667b9c50aac9b3ab6e042da229ae460b2380066fdf451fdbbad5`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 801–832. The root counter measured 1492098.345 µs, compared with 1,492,101 µs of draw time (1.78 ppm).

## Frame cadence

**Peak confirmed live render 47.186 ms** at frame 816 (GSReactionDiffusion). **0/945 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 24.668 ms. Live render averaged 38.932 ms.

Among complete live readouts, the largest window mean was 46.628 ms/frame at frames 801–832. The dominant scope averaged 39.343 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–946 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 947–960, with uncertain membership. Their largest render was 32.663 ms and 0 exceeded the budget. The live denominator lies between 945 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 47.186, spilled 🟢 0/945 (0.00%).

README setup frame 1: 24.668 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    24.352 ms 14.61 Mcyc 100.0%
  pov_preserve_half      135.24 us  81.1 kcyc   0.6% x0.97 140 us/call
  grd_render             24.064 ms 14.44 Mcyc  98.8% x1 24064 us/call
    grd_shader_draw       9.946 ms  5.97 Mcyc  40.8% x6 1658 us/call
    grd_cull_flags       711.28 us 426.8 kcyc   2.9% x6 119 us/call
    grd_publish           1.958 ms  1.18 Mcyc   8.0% x6 326 us/call
    grd_physics           5.145 ms  3.09 Mcyc  21.1% x6 858 us/call
    grd_pigment           3.067 ms  1.84 Mcyc  12.6% x1 3067 us/call
    rd_lattice_basis     130.79 us  78.5 kcyc   0.5% x1 131 us/call
    grd_color_noise       2.337 ms  1.40 Mcyc   9.6% x1 2337 us/call
  rd_timeline_step        52.25 us  31.4 kcyc   0.2% x1 52.3 us/call
  canvas_clear            88.56 us  53.1 kcyc   0.4% x1 88.6 us/call
  canvas_buffer_wait       0.13 us     76 cyc   0.0% x1 0.13 us/call
grd_color_palette        424.09 us 254.5 kcyc   1.7% x0.03 13571 us/call
grd_seed_reaction        429.32 us 257.6 kcyc   1.8% x0.03 13738 us/call
```

Wall min/mean/max was 19.625/24.351/26.994 ms. Render averaged 24.352 ms per frame, including 24.064 ms in `grd_render`.

This regime leaves 38.148 ms of mean render headroom per display window.

### Window containing the live peak (frames 801–832)

```
frame                    46.628 ms 27.98 Mcyc 100.0%
  pov_preserve_half      133.37 us  80.0 kcyc   0.3% x1 133 us/call
  grd_render             46.346 ms 27.81 Mcyc  99.4% x1 46346 us/call
    grd_shader_draw      27.604 ms 16.56 Mcyc  59.2% x6 4601 us/call
    grd_cull_flags       720.49 us 432.3 kcyc   1.5% x6 120 us/call
    grd_publish           4.681 ms  2.81 Mcyc  10.0% x6 780 us/call
    grd_physics           5.161 ms  3.10 Mcyc  11.1% x6 860 us/call
    grd_pigment           4.203 ms  2.52 Mcyc   9.0% x1 4203 us/call
    rd_lattice_basis     133.12 us  79.9 kcyc   0.3% x1 133 us/call
    grd_color_noise       2.344 ms  1.41 Mcyc   5.0% x1 2344 us/call
  rd_timeline_step        51.38 us  30.8 kcyc   0.1% x1 51.4 us/call
  canvas_clear            86.28 us  51.8 kcyc   0.2% x1 86.3 us/call
  canvas_buffer_wait       0.12 us     75 cyc   0.0% x1 0.12 us/call
```

Wall min/mean/max was 45.466/46.628/47.186 ms. Render averaged 46.628 ms per frame, including 46.346 ms in `grd_render`.

This regime leaves 15.872 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    27.830 ms 16.70 Mcyc 100.0%
  pov_preserve_half      136.54 us  81.9 kcyc   0.5% x1 137 us/call
  grd_render             27.546 ms 16.53 Mcyc  99.0% x1 27546 us/call
    grd_shader_draw      12.932 ms  7.76 Mcyc  46.5% x6 2155 us/call
    grd_cull_flags       715.58 us 429.3 kcyc   2.6% x6 119 us/call
    grd_publish           1.949 ms  1.17 Mcyc   7.0% x6 325 us/call
    grd_physics           5.145 ms  3.09 Mcyc  18.5% x6 857 us/call
    grd_pigment           3.570 ms  2.14 Mcyc  12.8% x1 3570 us/call
    rd_lattice_basis     129.95 us  78.0 kcyc   0.5% x1 130 us/call
    grd_color_noise       2.339 ms  1.40 Mcyc   8.4% x1 2339 us/call
  rd_timeline_step        50.99 us  30.6 kcyc   0.2% x1 51.0 us/call
  canvas_clear            84.70 us  50.8 kcyc   0.3% x1 84.7 us/call
  canvas_buffer_wait       0.12 us     75 cyc   0.0% x1 0.12 us/call
```

Wall min/mean/max was 24.089/27.830/32.046 ms. Render averaged 27.830 ms per frame, including 27.546 ms in `grd_render`.

This regime leaves 34.670 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 801–832, measured over the separate 1,999,763 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.9/frame  0.475/1.748/22.738 us  3.22%
  isr_pack           144.0/frame  5.923/7.445/20.998 us  1.71%
  isr_dma_submit     144.0/frame  0.626/0.921/5.678 us  0.21%
```

Inclusive wake time leaves 60.487 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `grd_render` took 46.346 ms/frame (99.4% of draw time).
2. `grd_shader_draw` took 27.604 ms/frame (59.2% of draw time).
3. `grd_physics` took 5.161 ms/frame (11.1% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the effect-local physics/shading functions and shared scan/filter kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=GSReactionDiffusion` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
