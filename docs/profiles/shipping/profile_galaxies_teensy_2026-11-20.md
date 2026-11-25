# Galaxies on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/galaxies_ship.log`. Supersedes `profile_galaxies_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Galaxies at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Galaxies profile 66 16` from the captured revision |

Captured 2026-11-20 16:13, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `079a0adbfe5dfd2e969b9958c82e807356c89df946373bb362f91c08d1239672`.

Profile image size `FLASH: code:92632, data:149444, headers:8804 free for files:1780736 / RAM1: variables:313952, code:23640, padding:9128 free for local variables:177568 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `c053fb917266d842f25d3b81bf1c2bb1cefd40d46cb9dc6d23a0c57d7c9be157`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 865–880. The root counter measured 687612.445 µs, compared with 687,613 µs of draw time (0.81 ppm).

## Frame cadence

**Peak confirmed live render 43.659 ms** at frame 879 (Galaxies). **0/949 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 1.470 ms. Live render averaged 24.752 ms.

Among complete live readouts, the largest window mean was 42.975 ms/frame at frames 865–880. The dominant scope averaged 15.827 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–950 are confirmed live. The raw file contains 1,032 frames, including 75 confirmed after the boundary.

Boundary resolution leaves 7 frame(s), 951–957, with uncertain membership. Their largest render was 43.183 ms and 0 exceeded the budget. The live denominator lies between 949 and 956. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 43.659, spilled 🟢 0/949 (0.00%).

README setup frame 1: 1.470 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                     1.522 ms 913.1 kcyc 100.0%
  pov_preserve_half      140.11 us  84.1 kcyc   9.2% x0.94 149 us/call
  gx_draw_particles      806.38 us 483.8 kcyc  53.0% x1 806 us/call
  gx_particle_step       395.08 us 237.0 kcyc  26.0% x1 395 us/call
  gx_timeline_step        56.15 us  33.7 kcyc   3.7% x1 56.2 us/call
  canvas_clear           101.21 us  60.7 kcyc   6.7% x1 101 us/call
  canvas_buffer_wait       1.23 us    738 cyc   0.1% x1 1.23 us/call
```

Wall min/mean/max was 1.164/1.521/1.846 ms. Render averaged 1.521 ms per frame, including 0.395 ms in `gx_particle_step`.

This regime leaves 60.979 ms of mean render headroom per display window.

### Window containing the live peak (frames 865–880)

```
frame                    42.976 ms 25.79 Mcyc 100.0%
  pov_preserve_half      143.87 us  86.3 kcyc   0.3% x1 144 us/call
  gx_draw_particles      15.775 ms  9.47 Mcyc  36.7% x1 15775 us/call
  gx_particle_step       26.875 ms 16.12 Mcyc  62.5% x1 26875 us/call
  gx_timeline_step        56.42 us  33.9 kcyc   0.1% x1 56.4 us/call
  canvas_clear            97.86 us  58.7 kcyc   0.2% x1 97.9 us/call
  canvas_buffer_wait       1.07 us    643 cyc   0.0% x1 1.07 us/call
```

Wall min/mean/max was 42.402/42.975/43.660 ms. Render averaged 42.975 ms per frame, including 26.875 ms in `gx_particle_step`.

This regime leaves 19.525 ms of mean render headroom per display window.

### Final complete in-show window (frames 929–944)

```
frame                    41.587 ms 24.95 Mcyc 100.0%
  pov_preserve_half      143.32 us  86.0 kcyc   0.3% x1 143 us/call
  gx_draw_particles      14.464 ms  8.68 Mcyc  34.8% x1 14464 us/call
  gx_particle_step       26.808 ms 16.08 Mcyc  64.5% x1 26808 us/call
  gx_timeline_step        53.91 us  32.3 kcyc   0.1% x1 53.9 us/call
  canvas_clear            94.79 us  56.9 kcyc   0.2% x1 94.8 us/call
  canvas_buffer_wait       1.07 us    644 cyc   0.0% x1 1.07 us/call
```

Wall min/mean/max was 40.900/41.587/42.556 ms. Render averaged 41.586 ms per frame, including 26.808 ms in `gx_particle_step`.

This regime leaves 20.914 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 865–880, measured over the separate 999,588 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.6/frame  0.480/1.595/17.986 us  2.94%
  isr_pack           143.9/frame  4.755/6.265/9.106 us  1.44%
  isr_dma_submit     143.9/frame  0.803/0.915/4.840 us  0.21%
```

Inclusive wake time leaves 60.663 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `gx_particle_step` took 26.875 ms/frame (62.5% of draw time).
2. `gx_draw_particles` took 15.775 ms/frame (36.7% of draw time).
3. `pov_preserve_half` took 0.144 ms/frame (0.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared scan-shapes, spherical-math and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Galaxies` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
