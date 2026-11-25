# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/bzreactiondiffusion_ship.log`. Supersedes `profile_bzreactiondiffusion_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | BZReactionDiffusion at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh BZReactionDiffusion profile 66 32` from the captured revision |

Captured 2026-11-20 15:54, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `1a39143b38c64551125a0e9539d39dab3129a8a6a917644d72de732a7d726608`.

Profile image size `FLASH: code:64904, data:156588, headers:8908 free for files:1801216 / RAM1: variables:314048, code:17304, padding:15464 free for local variables:177472 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `c12a686dd43ffcf94270da90da69be74e433cfbe4267312ca4e575573b25f5e8`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 65–96. The root counter measured 1601118.092 µs, compared with 1,601,119 µs of draw time (0.57 ppm).

## Frame cadence

**Peak confirmed live render 50.177 ms** at frame 604 (BZReactionDiffusion). **0/950 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 87.527 ms. Live render averaged 49.902 ms.

Among complete live readouts, the largest window mean was 50.034 ms/frame at frames 65–96. The dominant scope averaged 49.657 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–951 are confirmed live. The raw file contains 1,032 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 8 frame(s), 952–959, with uncertain membership. Their largest render was 50.093 ms and 0 exceeded the budget. The live denominator lies between 950 and 958. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 50.177, spilled 🟢 0/950 (0.00%).

README setup frame 1: 87.527 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    50.435 ms 30.26 Mcyc 100.0%
  pov_preserve_half      131.79 us  79.1 kcyc   0.3% x0.97 136 us/call
  bz_render              50.162 ms 30.10 Mcyc  99.5% x1 50162 us/call
    bz_raster            45.882 ms 27.53 Mcyc  91.0% x2 22941 us/call
    bz_physics            4.179 ms  2.51 Mcyc   8.3% x2 2089 us/call
    rd_lattice_basis      95.13 us  57.1 kcyc   0.2% x1 95.1 us/call
  rd_timeline_step        32.53 us  19.5 kcyc   0.1% x1 32.5 us/call
  canvas_clear            91.68 us  55.0 kcyc   0.2% x1 91.7 us/call
  canvas_buffer_wait       1.09 us    655 cyc   0.0% x1 1.09 us/call
```

Wall min/mean/max was 47.603/50.435/87.527 ms. Render averaged 50.434 ms per frame, including 50.162 ms in `bz_render`.

This regime leaves 12.066 ms of mean render headroom per display window.

### Window containing the live peak (frames 577–608)

```
frame                    49.990 ms 29.99 Mcyc 100.0%
  pov_preserve_half      135.57 us  81.3 kcyc   0.3% x1 136 us/call
  bz_render              49.722 ms 29.83 Mcyc  99.5% x1 49722 us/call
    bz_raster            45.481 ms 27.29 Mcyc  91.0% x2 22740 us/call
    bz_physics            4.142 ms  2.48 Mcyc   8.3% x2 2071 us/call
    rd_lattice_basis      94.62 us  56.8 kcyc   0.2% x1 94.6 us/call
  rd_timeline_step        31.18 us  18.7 kcyc   0.1% x1 31.2 us/call
  canvas_clear            87.55 us  52.5 kcyc   0.2% x1 87.5 us/call
  canvas_buffer_wait       1.09 us    652 cyc   0.0% x1 1.09 us/call
```

Wall min/mean/max was 49.824/49.989/50.178 ms. Render averaged 49.989 ms per frame, including 49.722 ms in `bz_render`.

This regime leaves 12.511 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    49.813 ms 29.89 Mcyc 100.0%
  pov_preserve_half      135.03 us  81.0 kcyc   0.3% x1 135 us/call
  bz_render              49.548 ms 29.73 Mcyc  99.5% x1 49548 us/call
    bz_raster            45.289 ms 27.17 Mcyc  90.9% x2 22645 us/call
    bz_physics            4.159 ms  2.50 Mcyc   8.3% x2 2080 us/call
    rd_lattice_basis      94.78 us  56.9 kcyc   0.2% x1 94.8 us/call
  rd_timeline_step        29.22 us  17.5 kcyc   0.1% x1 29.2 us/call
  canvas_clear            87.17 us  52.3 kcyc   0.2% x1 87.2 us/call
  canvas_buffer_wait       1.09 us    653 cyc   0.0% x1 1.09 us/call
```

Wall min/mean/max was 49.697/49.812/49.913 ms. Render averaged 49.811 ms per frame, including 49.548 ms in `bz_render`.

This regime leaves 12.689 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 577–608, measured over the separate 2,000,036 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.1/frame  0.476/1.760/13.871 us  3.24%
  isr_pack           144.0/frame  7.166/7.748/10.288 us  1.78%
  isr_dma_submit     144.0/frame  0.611/0.913/1.041 us  0.21%
```

Inclusive wake time leaves 60.475 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `bz_render` took 49.722 ms/frame (99.5% of draw time).
2. `bz_raster` took 45.481 ms/frame (91.0% of draw time).
3. `bz_physics` took 4.142 ms/frame (8.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the effect-local physics/shading functions and shared scan/filter kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=BZReactionDiffusion` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
