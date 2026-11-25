# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/gnomonicstars_ship.log`. Supersedes `profile_gnomonicstars_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | GnomonicStars at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh GnomonicStars profile 66 32` from the captured revision |

Captured 2026-11-20 15:57, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `55729c8e1afc9770ebafff7fbb38155d17acdc29fcaf914919eb0ebabdff798d`.

Profile image size `FLASH: code:75928, data:148964, headers:8580 free for files:1798144 / RAM1: variables:314048, code:20824, padding:11944 free for local variables:177472 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `bbd97291fc8ed1b84f584860bcfb998e0def972b494be3936ca958d05ad36fab`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 193–224. The root counter measured 592063.398 µs, compared with 592,064 µs of draw time (1.02 ppm).

## Frame cadence

**Peak confirmed live render 28.079 ms** at frame 787 (GnomonicStars). **0/928 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 15.235 ms. Live render averaged 12.011 ms.

Among complete live readouts, the largest window mean was 18.500 ms/frame at frames 193–224. The dominant scope averaged 11.733 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–929 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 31 frame(s), 930–960, with uncertain membership. Their largest render was 23.047 ms and 0 exceeded the budget. The live denominator lies between 928 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 28.079, spilled 🟢 0/928 (0.00%).

README setup frame 1: 15.235 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    10.485 ms  6.29 Mcyc 100.0%
  pov_preserve_half      137.59 us  82.6 kcyc   1.3% x0.97 142 us/call
  gn_draw_stars          10.123 ms  6.07 Mcyc  96.5% x1 10123 us/call
    gn_star_scan          8.481 ms  5.09 Mcyc  80.9% x252.81 33.5 us/call
  gn_spiral_build         17.50 us  10.5 kcyc   0.2% x0.03 560 us/call
  gn_timeline_step        62.91 us  37.7 kcyc   0.6% x1 62.9 us/call
    animation_mobius_step  17.45 us  10.5 kcyc   0.2% x1 17.4 us/call
  canvas_clear            98.42 us  59.1 kcyc   0.9% x1 98.4 us/call
  canvas_buffer_wait       2.12 us   1.3 kcyc   0.0% x1 2.12 us/call
```

Wall min/mean/max was 4.433/10.485/15.239 ms. Render averaged 10.483 ms per frame, including 10.123 ms in `gn_draw_stars`.

This regime leaves 52.017 ms of mean render headroom per display window.

### Window containing the live peak (frames 769–800)

```
frame                    14.918 ms  8.95 Mcyc 100.0%
  pov_preserve_half      139.46 us  83.7 kcyc   0.9% x1 139 us/call
  gn_draw_stars          14.576 ms  8.75 Mcyc  97.7% x1 14576 us/call
    gn_star_scan         12.990 ms  7.79 Mcyc  87.1% x241.78 53.7 us/call
  gn_timeline_step        61.34 us  36.8 kcyc   0.4% x1 61.3 us/call
    animation_mobius_step  17.95 us  10.8 kcyc   0.1% x1 18.0 us/call
  canvas_clear            93.57 us  56.1 kcyc   0.6% x1 93.6 us/call
  canvas_buffer_wait       2.21 us   1.3 kcyc   0.0% x1 2.21 us/call
```

Wall min/mean/max was 5.648/14.917/28.081 ms. Render averaged 14.916 ms per frame, including 14.576 ms in `gn_draw_stars`.

This regime leaves 47.584 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    12.635 ms  7.58 Mcyc 100.0%
  pov_preserve_half      143.14 us  85.9 kcyc   1.1% x1 143 us/call
  gn_draw_stars          12.289 ms  7.37 Mcyc  97.3% x1 12289 us/call
    gn_star_scan         10.704 ms  6.42 Mcyc  84.7% x233.34 45.9 us/call
  gn_timeline_step        55.32 us  33.2 kcyc   0.4% x1 55.3 us/call
    animation_mobius_step  16.38 us   9.8 kcyc   0.1% x1 16.4 us/call
  canvas_clear            95.35 us  57.2 kcyc   0.8% x1 95.3 us/call
  canvas_buffer_wait       1.92 us   1.1 kcyc   0.0% x1 1.92 us/call
```

Wall min/mean/max was 7.708/12.634/17.197 ms. Render averaged 12.633 ms per frame, including 12.289 ms in `gn_draw_stars`.

This regime leaves 49.867 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 769–800, measured over the separate 1,992,673 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1147.8/frame  0.501/1.488/18.385 us  2.74%
  isr_pack           143.5/frame  4.536/5.643/9.223 us  1.30%
  isr_dma_submit     143.5/frame  0.698/0.906/0.973 us  0.20%
```

Inclusive wake time leaves 60.788 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `gn_draw_stars` took 14.576 ms/frame (97.7% of draw time).
2. `gn_star_scan` took 12.990 ms/frame (87.1% of draw time).
3. `pov_preserve_half` took 0.139 ms/frame (0.9% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared scan-shapes, spherical-math and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=GnomonicStars` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
