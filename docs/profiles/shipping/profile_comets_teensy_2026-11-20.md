# Comets on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 120-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/comets_ship.log`. Supersedes `profile_comets_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Comets at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 120 s slot, 126 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Comets profile 126 16` from the captured revision |

Captured 2026-11-20 15:19, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `ed6345dc687eaa0d334575908e48a1d1f7e6513900f5098a7da1cdc11cd34ab7`.

Profile image size `FLASH: code:86976, data:148416, headers:8320 free for files:1787904 / RAM1: variables:313952, code:26408, padding:6360 free for local variables:177568 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `beb0611c16c6f99fe2fde29acb211b018b2bdd8547f55e7b435c25374129c725`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 513–528. The root counter measured 360540.135 µs, compared with 360,539 µs of draw time (3.15 ppm).

## Frame cadence

**Peak confirmed live render 29.874 ms** at frame 647 (5). **0/1,903 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 0.879 ms. Live render averaged 15.043 ms.

Among complete live readouts, the largest window mean was 22.532 ms/frame at frames 513–528. The dominant scope averaged 14.232 ms/frame across those readouts.

The shipping ELF assigns 960 revolutions (120 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 120 s. Frames 2–1904 are confirmed live. The raw file contains 1,993 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 1905–1920, with uncertain membership. Their largest render was 16.013 ms and 0 exceeded the budget. The live denominator lies between 1,903 and 1,919. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 29.874 (12), spilled 🟢 0/1903 (0.00%).

README setup frame 1: 0.879 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                     2.114 ms  1.27 Mcyc 100.0%
  pov_preserve_half      140.49 us  84.3 kcyc   6.6% x0.94 150 us/call
  cm_draw_trail           1.772 ms  1.06 Mcyc  83.8% x1 1772 us/call
  cm_wipe_rebake           0.09 us     56 cyc   0.0% x1 0.09 us/call
  cm_timeline_step        96.53 us  57.9 kcyc   4.6% x1 96.5 us/call
  canvas_clear            94.24 us  56.5 kcyc   4.5% x1 94.2 us/call
  canvas_buffer_wait       2.81 us   1.7 kcyc   0.1% x1 2.81 us/call
```

Wall min/mean/max was 0.882/2.114/3.039 ms. Render averaged 2.111 ms per frame, including 1.772 ms in `cm_draw_trail`.

This regime leaves 60.389 ms of mean render headroom per display window.

### Window containing the live peak (frames 641–656)

```
frame                    21.952 ms 13.17 Mcyc 100.0%
  pov_preserve_half      136.45 us  81.9 kcyc   0.6% x1 136 us/call
  cm_draw_trail          19.626 ms 11.78 Mcyc  89.4% x1 19626 us/call
  cm_wipe_rebake          1.925 ms  1.15 Mcyc   8.8% x1 1925 us/call
  cm_timeline_step       169.90 us 101.9 kcyc   0.8% x1 170 us/call
  canvas_clear            86.24 us  51.7 kcyc   0.4% x1 86.2 us/call
  canvas_buffer_wait       1.60 us    962 cyc   0.0% x1 1.60 us/call
```

Wall min/mean/max was 12.418/21.951/29.876 ms. Render averaged 21.950 ms per frame, including 19.626 ms in `cm_draw_trail`.

This regime leaves 40.550 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 513–528)

```
frame                    22.534 ms 13.52 Mcyc 100.0%
  pov_preserve_half      138.54 us  83.1 kcyc   0.6% x1 139 us/call
  cm_draw_trail          20.178 ms 12.11 Mcyc  89.5% x1 20178 us/call
  cm_wipe_rebake          1.949 ms  1.17 Mcyc   8.7% x1 1949 us/call
  cm_timeline_step       174.15 us 104.5 kcyc   0.8% x1 174 us/call
  canvas_clear            85.75 us  51.4 kcyc   0.4% x1 85.7 us/call
  canvas_buffer_wait       1.55 us    929 cyc   0.0% x1 1.55 us/call
```

Wall min/mean/max was 20.559/22.533/24.844 ms. Render averaged 22.532 ms per frame, including 20.178 ms in `cm_draw_trail`.

This regime leaves 39.968 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 625–640)

```
frame                    19.976 ms 11.99 Mcyc 100.0%
  pov_preserve_half      139.07 us  83.4 kcyc   0.7% x1 139 us/call
  cm_draw_trail          19.601 ms 11.76 Mcyc  98.1% x1 19601 us/call
  cm_wipe_rebake           0.09 us     54 cyc   0.0% x1 0.09 us/call
  cm_timeline_step       136.79 us  82.1 kcyc   0.7% x1 137 us/call
  canvas_clear            85.87 us  51.5 kcyc   0.4% x1 85.9 us/call
  canvas_buffer_wait       1.67 us   1.0 kcyc   0.0% x1 1.67 us/call
```

Wall min/mean/max was 18.003/19.975/23.438 ms. Render averaged 19.974 ms per frame, including 19.601 ms in `cm_draw_trail`.

This regime leaves 42.526 ms of mean render headroom per display window.

### Final complete in-show window (frames 1889–1904)

```
frame                    15.746 ms  9.45 Mcyc 100.0%
  pov_preserve_half      141.97 us  85.2 kcyc   0.9% x1 142 us/call
  cm_draw_trail          15.380 ms  9.23 Mcyc  97.7% x1 15380 us/call
  cm_wipe_rebake           0.09 us     55 cyc   0.0% x1 0.09 us/call
  cm_timeline_step       129.14 us  77.5 kcyc   0.8% x1 129 us/call
  canvas_clear            86.02 us  51.6 kcyc   0.5% x1 86.0 us/call
  canvas_buffer_wait       1.67 us   1.0 kcyc   0.0% x1 1.67 us/call
```

Wall min/mean/max was 12.849/15.746/18.093 ms. Render averaged 15.744 ms per frame, including 15.380 ms in `cm_draw_trail`.

This regime leaves 46.756 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| 4 | held preset | 9 | 513–528 | — | — | 20.178 | 22.532 | 25.910 | 16 | 0/160 |
| 10 | held preset | 9 | 1489–1504 | — | — | 19.820 | 20.188 | 26.610 | 16 | 0/160 |
| 12 | held preset | 8 | 1857–1872 | — | — | 19.721 | 20.084 | 22.914 | 16 | 0/145 |
| 5 | held preset | 9 | 641–656 | — | — | 19.626 | 21.950 | 29.874 | 16 | 0/160 |
| 2 | held preset | 9 | 289–304 | — | — | 19.346 | 19.703 | 24.750 | 16 | 0/160 |
| 8 | held preset | 9 | 1233–1248 | — | — | 19.115 | 19.475 | 22.054 | 16 | 0/160 |
| 9 | held preset | 9 | 1313–1328 | — | — | 18.637 | 21.015 | 25.293 | 16 | 0/160 |
| 11 | held preset | 9 | 1633–1648 | — | — | 17.690 | 20.030 | 26.932 | 16 | 0/160 |
| 6 | held preset | 9 | 801–816 | — | — | 14.698 | 16.978 | 20.574 | 16 | 0/160 |
| 3 | held preset | 9 | 433–448 | — | — | 14.663 | 15.011 | 18.787 | 16 | 0/160 |
| 7 | held preset | 9 | 1089–1104 | — | — | 12.339 | 12.675 | 15.467 | 16 | 0/160 |
| 1 | held preset | 8 | 113–128 | — | — | 12.077 | 12.399 | 14.571 | 16 | 0/158 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 641–656, measured over the separate 994,412 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1145.6/frame  0.533/1.506/16.491 us  2.77%
  isr_pack           143.2/frame  4.536/5.570/8.968 us  1.28%
  isr_dma_submit     143.2/frame  0.706/0.910/4.198 us  0.20%
```

Inclusive wake time leaves 60.769 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `cm_draw_trail` took 19.626 ms/frame (89.4% of draw time).
2. `cm_wipe_rebake` took 1.925 ms/frame (8.8% of draw time).
3. `cm_timeline_step` took 0.170 ms/frame (0.8% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared Plot raster, column-walk and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Comets` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
