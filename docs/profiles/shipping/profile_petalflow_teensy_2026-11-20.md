# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 30-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/petalflow_ship.log`. Supersedes `profile_petalflow_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | PetalFlow at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 30 s slot, 36 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh PetalFlow profile 36 32` from the captured revision |

Captured 2026-11-20 16:21, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `4bd9c8b481e18e36c8deef478c7e94fd55de0fd9505079779707506a9c8ab9cb`.

Profile image size `FLASH: code:113384, data:145768, headers:9136 free for files:1763328 / RAM1: variables:312992, code:29720, padding:3048 free for local variables:178528 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `176dbf7278bd17e0d099981c0b59c9161f297b162006683977f601e1d52d117a`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 321–352. The root counter measured 748669.983 µs, compared with 748,671 µs of draw time (1.36 ppm).

## Frame cadence

**Peak confirmed live render 25.346 ms** at frame 307 (PetalFlow). **0/459 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 19.772 ms. Live render averaged 20.915 ms.

Among complete live readouts, the largest window mean was 23.395 ms/frame at frames 321–352. The dominant scope averaged 20.776 ms/frame across those readouts.

The shipping ELF assigns 240 revolutions (30 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 30 s. Frames 2–460 are confirmed live. The raw file contains 553 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 20 frame(s), 461–480, with uncertain membership. Their largest render was 21.569 ms and 0 exceeded the budget. The live denominator lies between 459 and 479. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 25.346, spilled 🟢 0/459 (0.00%).

README setup frame 1: 19.772 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    19.131 ms 11.48 Mcyc 100.0%
  pov_preserve_half      134.35 us  80.6 kcyc   0.7% x0.97 139 us/call
  pf_draw_rings          18.834 ms 11.30 Mcyc  98.4% x1 18834 us/call
    pf_ring_scan         17.883 ms 10.73 Mcyc  93.5% x23.34 766 us/call
    pf_ring_build        850.56 us 510.3 kcyc   4.4% x23.34 36.4 us/call
  pf_timeline_step        45.02 us  27.0 kcyc   0.2% x1 45.0 us/call
  canvas_clear            91.87 us  55.1 kcyc   0.5% x1 91.9 us/call
  canvas_buffer_wait       1.17 us    699 cyc   0.0% x1 1.17 us/call
```

Wall min/mean/max was 17.088/19.131/21.775 ms. Render averaged 19.130 ms per frame, including 18.834 ms in `pf_draw_rings`.

This regime leaves 43.370 ms of mean render headroom per display window.

### Window containing the live peak (frames 289–320)

```
frame                    23.347 ms 14.01 Mcyc 100.0%
  pov_preserve_half      135.66 us  81.4 kcyc   0.6% x1 136 us/call
  pf_draw_rings          23.053 ms 13.83 Mcyc  98.7% x1 23053 us/call
    pf_ring_scan         22.110 ms 13.27 Mcyc  94.7% x23.34 947 us/call
    pf_ring_build        846.31 us 507.8 kcyc   3.6% x23.34 36.3 us/call
  pf_timeline_step        44.54 us  26.7 kcyc   0.2% x1 44.5 us/call
  canvas_clear            88.14 us  52.9 kcyc   0.4% x1 88.1 us/call
  canvas_buffer_wait       1.30 us    778 cyc   0.0% x1 1.30 us/call
```

Wall min/mean/max was 20.787/23.346/25.348 ms. Render averaged 23.345 ms per frame, including 23.053 ms in `pf_draw_rings`.

This regime leaves 39.155 ms of mean render headroom per display window.

### Final complete in-show window (frames 417–448)

```
frame                    20.637 ms 12.38 Mcyc 100.0%
  pov_preserve_half      141.12 us  84.7 kcyc   0.7% x1 141 us/call
  pf_draw_rings          20.340 ms 12.20 Mcyc  98.6% x1 20340 us/call
    pf_ring_scan         19.393 ms 11.64 Mcyc  94.0% x23.34 831 us/call
    pf_ring_build        843.43 us 506.1 kcyc   4.1% x23.34 36.1 us/call
  pf_timeline_step        44.52 us  26.7 kcyc   0.2% x1 44.5 us/call
  canvas_clear            87.85 us  52.7 kcyc   0.4% x1 87.9 us/call
  canvas_buffer_wait       1.46 us    878 cyc   0.0% x1 1.46 us/call
```

Wall min/mean/max was 19.399/20.637/22.221 ms. Render averaged 20.636 ms per frame, including 20.340 ms in `pf_draw_rings`.

This regime leaves 41.864 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 289–320, measured over the separate 2,003,427 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1154.0/frame  0.546/1.478/26.860 us  2.72%
  isr_pack           144.2/frame  4.621/5.468/8.441 us  1.25%
  isr_dma_submit     144.2/frame  0.695/0.918/9.230 us  0.21%
```

Inclusive wake time leaves 60.800 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `pf_draw_rings` took 23.053 ms/frame (98.7% of draw time).
2. `pf_ring_scan` took 22.110 ms/frame (94.7% of draw time).
3. `pf_ring_build` took 0.846 ms/frame (3.6% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared scan-shapes, spherical-math and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=PetalFlow` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
