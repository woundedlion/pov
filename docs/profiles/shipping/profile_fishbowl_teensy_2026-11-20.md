# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/fishbowl_ship.log`. Supersedes `profile_fishbowl_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Fishbowl at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Fishbowl profile 66 32` from the captured revision |

Captured 2026-11-20 16:12, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `974f2b224879ecc38a75eeb19101d1dfc201d868e4871e79ece449bf7ee355d9`.

Profile image size `FLASH: code:139760, data:148484, headers:8716 free for files:1734656 / RAM1: variables:314048, code:37368, padding:28168 free for local variables:144704 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `ab96d7d71c05ec4df9abd605c58b3f6f834309af4181463f4b53621588cb8cff`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 481–512. The root counter measured 724805.400 µs, compared with 724,804 µs of draw time (1.93 ppm).

## Frame cadence

**Peak confirmed live render 27.892 ms** at frame 215 (Fishbowl). **0/941 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 1.653 ms. Live render averaged 21.518 ms.

Among complete live readouts, the largest window mean was 24.245 ms/frame at frames 193–224. The dominant scope averaged 12.223 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–942 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 18 frame(s), 943–960, with uncertain membership. Their largest render was 25.213 ms and 0 exceeded the budget. The live denominator lies between 941 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 27.892, spilled 🟢 0/941 (0.00%).

README setup frame 1: 1.653 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                     4.382 ms  2.63 Mcyc 100.0%
  pov_preserve_half      140.01 us  84.0 kcyc   3.2% x0.97 145 us/call
  fish_multiline_draw     1.873 ms  1.12 Mcyc  42.7% x1 1873 us/call
  fish_build_vertices     1.872 ms  1.12 Mcyc  42.7% x1 1872 us/call
  fish_noise_prepare       4.87 us   2.9 kcyc   0.1% x1 4.87 us/call
  fish_timeline_step     222.62 us 133.6 kcyc   5.1% x1 223 us/call
    animation_noise_clock   0.18 us    106 cyc   0.0% x1 0.18 us/call
  canvas_clear            91.97 us  55.2 kcyc   2.1% x1 92.0 us/call
  canvas_buffer_wait       1.99 us   1.2 kcyc   0.0% x1 1.99 us/call
```

Wall min/mean/max was 1.420/4.382/8.104 ms. Render averaged 4.380 ms per frame, including 1.872 ms in `fish_build_vertices`.

This regime leaves 58.120 ms of mean render headroom per display window.

### Window containing the live peak (frames 193–224)

```
frame                    24.247 ms 14.55 Mcyc 100.0%
  pov_preserve_half      138.90 us  83.3 kcyc   0.6% x1 139 us/call
  fish_multiline_draw    10.523 ms  6.31 Mcyc  43.4% x1 10523 us/call
  fish_build_vertices    12.611 ms  7.57 Mcyc  52.0% x1 12611 us/call
  fish_noise_prepare       5.50 us   3.3 kcyc   0.0% x1 5.50 us/call
  fish_timeline_step     230.21 us 138.1 kcyc   0.9% x1 230 us/call
    animation_noise_clock   0.18 us    106 cyc   0.0% x1 0.18 us/call
  canvas_clear            88.30 us  53.0 kcyc   0.4% x1 88.3 us/call
  canvas_buffer_wait       1.79 us   1.1 kcyc   0.0% x1 1.79 us/call
```

Wall min/mean/max was 21.672/24.246/27.893 ms. Render averaged 24.245 ms per frame, including 12.611 ms in `fish_build_vertices`.

This regime leaves 38.255 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    23.380 ms 14.03 Mcyc 100.0%
  pov_preserve_half      138.82 us  83.3 kcyc   0.6% x1 139 us/call
  fish_multiline_draw     9.648 ms  5.79 Mcyc  41.3% x1 9648 us/call
  fish_build_vertices    12.616 ms  7.57 Mcyc  54.0% x1 12616 us/call
  fish_noise_prepare       5.61 us   3.4 kcyc   0.0% x1 5.61 us/call
  fish_timeline_step     229.47 us 137.7 kcyc   1.0% x1 229 us/call
    animation_noise_clock   0.18 us    106 cyc   0.0% x1 0.18 us/call
  canvas_clear            88.14 us  52.9 kcyc   0.4% x1 88.1 us/call
  canvas_buffer_wait       2.24 us   1.3 kcyc   0.0% x1 2.24 us/call
```

Wall min/mean/max was 21.586/23.380/25.860 ms. Render averaged 23.378 ms per frame, including 12.616 ms in `fish_build_vertices`.

This regime leaves 39.122 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 193–224, measured over the separate 2,001,037 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.7/frame  0.550/1.485/24.003 us  2.73%
  isr_pack           144.1/frame  4.565/5.265/7.716 us  1.21%
  isr_dma_submit     144.1/frame  0.656/0.921/11.160 us  0.21%
```

Inclusive wake time leaves 60.794 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fish_build_vertices` took 12.611 ms/frame (52.0% of draw time).
2. `fish_multiline_draw` took 10.523 ms/frame (43.4% of draw time).
3. `fish_timeline_step` took 0.230 ms/frame (0.9% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared scan-shapes, spherical-math and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Fishbowl` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
