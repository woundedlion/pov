# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/ringspin_ship.log`. Supersedes `profile_ringspin_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | RingSpin at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RingSpin profile 66 32` from the captured revision |

Captured 2026-11-20 16:20, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `1eaaa480b95ca74d42b85625f27cac039428210618fc81846a37149a4c02eee3`.

Profile image size `FLASH: code:86712, data:148028, headers:8972 free for files:1787904 / RAM1: variables:313952, code:27336, padding:5432 free for local variables:177568 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `6cc8098c0b2c8fddb2e08385610bd6393ab6828347d481e3972e3c02c27f3b6b`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 737–768. The root counter measured 1259734.150 µs, compared with 1,259,734 µs of draw time (0.12 ppm).

## Frame cadence

**Peak confirmed live render 50.813 ms** at frame 761 (RingSpin). **0/948 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 3.257 ms. Live render averaged 36.473 ms.

Among complete live readouts, the largest window mean was 39.367 ms/frame at frames 737–768. The dominant scope averaged 36.563 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–949 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 11 frame(s), 950–960, with uncertain membership. Their largest render was 40.498 ms and 0 exceeded the budget. The live denominator lies between 948 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 50.813, spilled 🟢 0/948 (0.00%).

README setup frame 1: 3.257 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    23.996 ms 14.40 Mcyc 100.0%
  pov_preserve_half      133.80 us  80.3 kcyc   0.6% x0.97 138 us/call
  rs_draw_rings          23.699 ms 14.22 Mcyc  98.8% x1 23699 us/call
    rs_ring_scan         22.991 ms 13.79 Mcyc  95.8% x54.62 421 us/call
  rs_timeline_step        66.88 us  40.1 kcyc   0.3% x1 66.9 us/call
  canvas_clear            90.19 us  54.1 kcyc   0.4% x1 90.2 us/call
  canvas_buffer_wait       0.18 us    110 cyc   0.0% x1 0.18 us/call
```

Wall min/mean/max was 3.230/23.996/45.934 ms. Render averaged 23.996 ms per frame, including 23.699 ms in `rs_draw_rings`.

This regime leaves 38.504 ms of mean render headroom per display window.

### Window containing the live peak (frames 737–768)

```
frame                    39.367 ms 23.62 Mcyc 100.0%
  pov_preserve_half      140.65 us  84.4 kcyc   0.4% x1 141 us/call
  rs_draw_rings          39.070 ms 23.44 Mcyc  99.2% x1 39070 us/call
    rs_ring_scan         38.025 ms 22.81 Mcyc  96.6% x76 500 us/call
  rs_timeline_step        64.51 us  38.7 kcyc   0.2% x1 64.5 us/call
  canvas_clear            86.38 us  51.8 kcyc   0.2% x1 86.4 us/call
  canvas_buffer_wait       0.12 us     75 cyc   0.0% x1 0.12 us/call
```

Wall min/mean/max was 29.016/39.366/50.813 ms. Render averaged 39.367 ms per frame, including 39.070 ms in `rs_draw_rings`.

This regime leaves 23.133 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    35.751 ms 21.45 Mcyc 100.0%
  pov_preserve_half      140.84 us  84.5 kcyc   0.4% x1 141 us/call
  rs_draw_rings          35.454 ms 21.27 Mcyc  99.2% x1 35454 us/call
    rs_ring_scan         34.446 ms 20.67 Mcyc  96.4% x76 453 us/call
  rs_timeline_step        64.21 us  38.5 kcyc   0.2% x1 64.2 us/call
  canvas_clear            86.17 us  51.7 kcyc   0.2% x1 86.2 us/call
  canvas_buffer_wait       0.12 us     75 cyc   0.0% x1 0.12 us/call
```

Wall min/mean/max was 29.739/35.750/40.933 ms. Render averaged 35.751 ms per frame, including 35.454 ms in `rs_draw_rings`.

This regime leaves 26.749 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 737–768, measured over the separate 2,004,761 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1154.8/frame  0.445/1.658/16.340 us  3.05%
  isr_pack           144.3/frame  4.638/6.901/10.050 us  1.59%
  isr_dma_submit     144.3/frame  0.700/0.910/1.816 us  0.20%
```

Inclusive wake time leaves 60.594 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `rs_draw_rings` took 39.070 ms/frame (99.2% of draw time).
2. `rs_ring_scan` took 38.025 ms/frame (96.6% of draw time).
3. `pov_preserve_half` took 0.141 ms/frame (0.4% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared Plot raster, column-walk and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=RingSpin` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
