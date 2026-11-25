# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 105-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/riverstones_ship.log`. Supersedes `profile_riverstones_teensy_2026-11-18.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | RiverStones at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 105 s slot, 111 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh RiverStones profile 111 16` from the captured revision |

Captured 2026-11-20 15:27, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `6c3b888fa4cb46861ad336bf4fc439a24fdf18ac9f66ea124082e5d963e0cc2d`.

Profile image size `FLASH: code:125728, data:156732, headers:8356 free for files:1740800 / RAM1: variables:314080, code:30120, padding:2648 free for local variables:177440 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `2381a688aee20ff6d5e927a0b629735641cddf0358cfc640c16deec51b53ad22`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 993–1008. The root counter measured 914001.820 µs, compared with 914,003 µs of draw time (1.29 ppm).

## Frame cadence

**Peak confirmed live render 58.139 ms** at frame 1012 (river-stones → river-stones-2). **0/1,670 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 69.048 ms. Live render averaged 50.913 ms.

Among complete live readouts, the largest window mean was 57.124 ms/frame at frames 993–1008. The dominant scope averaged 45.069 ms/frame across those readouts.

The shipping ELF assigns 840 revolutions (105 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 105 s. Frames 2–1671 are confirmed live. The raw file contains 1,752 frames, including 77 confirmed after the boundary.

Boundary resolution leaves 4 frame(s), 1672–1675, with uncertain membership. Their largest render was 49.557 ms and 0 exceeded the budget. The live denominator lies between 1,670 and 1,674. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 58.139 (2), spilled 🟢 0/1670 (0.00%).

README setup frame 1: 69.048 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    44.286 ms 26.57 Mcyc 100.0%
  pov_preserve_half      131.99 us  79.2 kcyc   0.3% x0.94 141 us/call
  fx_shader_draw         38.685 ms 23.21 Mcyc  87.4% x1 38685 us/call
    fx_prepare_column     1.486 ms 891.4 kcyc   3.4% x144 10.3 us/call
  fx_prepare_frame        5.026 ms  3.02 Mcyc  11.3% x1 5026 us/call
  fx_advance              55.84 us  33.5 kcyc   0.1% x1 55.8 us/call
  fx_timeline_step        93.07 us  55.8 kcyc   0.2% x1 93.1 us/call
  canvas_clear            96.00 us  57.6 kcyc   0.2% x1 96.0 us/call
  canvas_buffer_wait       1.25 us    749 cyc   0.0% x1 1.25 us/call
```

Wall min/mean/max was 41.817/44.285/69.049 ms. Render averaged 44.285 ms per frame, including 38.685 ms in `fx_shader_draw`.

This regime leaves 18.215 ms of mean render headroom per display window.

### Window containing the live peak (frames 1009–1024)

```
frame                    56.890 ms 34.13 Mcyc 100.0%
  pov_preserve_half      135.40 us  81.2 kcyc   0.2% x1 135 us/call
  fx_shader_draw         50.797 ms 30.48 Mcyc  89.3% x1 50797 us/call
    fx_prepare_column     2.518 ms  1.51 Mcyc   4.4% x144 17.5 us/call
  fx_prepare_frame        5.329 ms  3.20 Mcyc   9.4% x1 5329 us/call
  fx_advance              64.28 us  38.6 kcyc   0.1% x1 64.3 us/call
  fx_timeline_step       273.77 us 164.3 kcyc   0.5% x1 274 us/call
  canvas_clear            89.94 us  54.0 kcyc   0.2% x1 89.9 us/call
  canvas_buffer_wait       1.38 us    828 cyc   0.0% x1 1.38 us/call
```

Wall min/mean/max was 54.905/56.890/58.140 ms. Render averaged 56.889 ms per frame, including 50.797 ms in `fx_shader_draw`.

This regime leaves 5.611 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 1153–1168)

```
frame                    54.234 ms 32.54 Mcyc 100.0%
  pov_preserve_half      135.61 us  81.4 kcyc   0.3% x1 136 us/call
  fx_shader_draw         48.346 ms 29.01 Mcyc  89.1% x1 48346 us/call
    fx_prepare_column     1.631 ms 978.4 kcyc   3.0% x144 11.3 us/call
  fx_prepare_frame        5.334 ms  3.20 Mcyc   9.8% x1 5334 us/call
  fx_advance              57.66 us  34.6 kcyc   0.1% x1 57.7 us/call
  fx_timeline_step        83.74 us  50.2 kcyc   0.2% x1 83.7 us/call
  canvas_clear            88.39 us  53.0 kcyc   0.2% x1 88.4 us/call
  canvas_buffer_wait       1.25 us    748 cyc   0.0% x1 1.25 us/call
```

Wall min/mean/max was 52.878/54.233/57.247 ms. Render averaged 54.233 ms per frame, including 48.346 ms in `fx_shader_draw`.

This regime leaves 8.267 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 993–1008)

```
frame                    57.125 ms 34.28 Mcyc 100.0%
  pov_preserve_half      138.68 us  83.2 kcyc   0.2% x1 139 us/call
  fx_shader_draw         51.040 ms 30.62 Mcyc  89.3% x1 51040 us/call
    fx_prepare_column     2.513 ms  1.51 Mcyc   4.4% x144 17.5 us/call
  fx_prepare_frame        5.324 ms  3.19 Mcyc   9.3% x1 5324 us/call
  fx_advance              64.09 us  38.5 kcyc   0.1% x1 64.1 us/call
  fx_timeline_step       270.96 us 162.6 kcyc   0.5% x1 271 us/call
  canvas_clear            89.44 us  53.7 kcyc   0.2% x1 89.4 us/call
  canvas_buffer_wait       1.25 us    748 cyc   0.0% x1 1.25 us/call
```

Wall min/mean/max was 56.082/57.125/57.827 ms. Render averaged 57.124 ms per frame, including 51.040 ms in `fx_shader_draw`.

This regime leaves 5.376 ms of mean render headroom per display window.

### Final complete in-show window (frames 1649–1664)

```
frame                    47.987 ms 28.79 Mcyc 100.0%
  pov_preserve_half      137.61 us  82.6 kcyc   0.3% x1 138 us/call
  fx_shader_draw         42.121 ms 25.27 Mcyc  87.8% x1 42121 us/call
    fx_prepare_column     1.617 ms 970.0 kcyc   3.4% x144 11.2 us/call
  fx_prepare_frame        5.315 ms  3.19 Mcyc  11.1% x1 5315 us/call
  fx_advance              56.28 us  33.8 kcyc   0.1% x1 56.3 us/call
  fx_timeline_step        78.93 us  47.4 kcyc   0.2% x1 78.9 us/call
  canvas_clear            88.19 us  52.9 kcyc   0.2% x1 88.2 us/call
  canvas_buffer_wait       1.31 us    788 cyc   0.0% x1 1.31 us/call
```

Wall min/mean/max was 47.486/47.986/49.202 ms. Render averaged 47.985 ms per frame, including 42.121 ms in `fx_shader_draw`.

This regime leaves 14.515 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| river-stones-2 | held preset | 36 | 1153–1168 | — | — | 48.346 | 54.233 | 58.139 | 16 | 0/1072 |
| river-stones | held preset | 36 | 193–208 | — | — | 48.154 | 54.047 | 54.582 | 16 | 0/598 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 1009–1024, measured over the separate 999,512 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.5/frame  0.491/1.801/18.445 us  3.32%
  isr_pack           143.9/frame  6.616/7.611/10.323 us  1.75%
  isr_dma_submit     143.9/frame  0.745/0.918/5.676 us  0.21%
```

Inclusive wake time leaves 60.425 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 50.797 ms/frame (89.3% of draw time).
2. `fx_prepare_frame` took 5.329 ms/frame (9.4% of draw time).
3. `fx_prepare_column` took 2.518 ms/frame (4.4% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=RiverStones` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
