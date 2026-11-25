# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 200-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/dreamballs_ship.log`. Supersedes `profile_dreamballs_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | DreamBalls at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 200 s slot, 206 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh DreamBalls profile 206 16` from the captured revision |

Captured 2026-11-20 15:36, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `bee7a3205e34307679ac737b377d90f568a2f87e16ef5ba1cfcada3ac05642fd`.

Profile image size `FLASH: code:153392, data:189376, headers:8464 free for files:1680384 / RAM1: variables:314144, code:35928, padding:29608 free for local variables:144608 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `edde953fcfd5f0c020014b807b772c457cadd1e429accd01c276bbc43a6b895c`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 2769–2784. The root counter measured 776045.560 µs, compared with 776,044 µs of draw time (2.01 ppm).

## Frame cadence

**Peak confirmed live render 51.441 ms** at frame 2752 (9). **0/3,187 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 51.219 ms. Live render averaged 28.275 ms.

Among complete live readouts, the largest window mean was 48.502 ms/frame at frames 2769–2784. The dominant scope averaged 27.969 ms/frame across those readouts.

The shipping ELF assigns 1,600 revolutions (200 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 200 s. Frames 2–3188 are confirmed live. The raw file contains 3,272 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 12 frame(s), 3189–3200, with uncertain membership. Their largest render was 34.469 ms and 0 exceeded the budget. The live denominator lies between 3,187 and 3,199. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 51.441 (10), spilled 🟢 0/3187 (0.00%).

README setup frame 1: 51.219 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    38.951 ms 23.37 Mcyc 100.0%
  pov_preserve_half      136.10 us  81.7 kcyc   0.3% x0.94 145 us/call
  db_timeline_step       38.697 ms 23.22 Mcyc  99.3% x1 38697 us/call
    db_draw              38.614 ms 23.17 Mcyc  99.1% x1 38614 us/call
      db_draw_scene      38.600 ms 23.16 Mcyc  99.1% x1 38600 us/call
        db_mesh_plot     38.050 ms 22.83 Mcyc  97.7% x18 2114 us/call
        db_displace      479.59 us 287.8 kcyc   1.2% x18 26.6 us/call
  canvas_clear            98.94 us  59.4 kcyc   0.3% x1 98.9 us/call
  canvas_buffer_wait       1.14 us    682 cyc   0.0% x1 1.14 us/call
```

Wall min/mean/max was 36.284/38.951/51.221 ms. Render averaged 38.950 ms per frame, including 38.697 ms in `db_timeline_step`.

This regime leaves 23.550 ms of mean render headroom per display window.

### Window containing the live peak (frames 2737–2752)

```
frame                    48.133 ms 28.88 Mcyc 100.0%
  pov_preserve_half      140.15 us  84.1 kcyc   0.3% x1 140 us/call
  db_timeline_step       47.883 ms 28.73 Mcyc  99.5% x1 47883 us/call
    db_draw              47.737 ms 28.64 Mcyc  99.2% x1 47737 us/call
      db_draw_scene      47.722 ms 28.63 Mcyc  99.1% x1 47722 us/call
        db_mesh_plot     46.348 ms 27.81 Mcyc  96.3% x6 7725 us/call
        db_displace       1.124 ms 674.6 kcyc   2.3% x6 187 us/call
  canvas_clear            91.68 us  55.0 kcyc   0.2% x1 91.7 us/call
  canvas_buffer_wait       0.83 us    501 cyc   0.0% x1 0.83 us/call
```

Wall min/mean/max was 45.492/48.132/51.442 ms. Render averaged 48.132 ms per frame, including 47.883 ms in `db_timeline_step`.

This regime leaves 14.368 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 2769–2784)

```
frame                    48.503 ms 29.10 Mcyc 100.0%
  pov_preserve_half      142.45 us  85.5 kcyc   0.3% x1 142 us/call
  db_timeline_step       48.250 ms 28.95 Mcyc  99.5% x1 48250 us/call
    db_draw              48.102 ms 28.86 Mcyc  99.2% x1 48102 us/call
      db_draw_scene      48.088 ms 28.85 Mcyc  99.1% x1 48088 us/call
        db_mesh_plot     46.707 ms 28.02 Mcyc  96.3% x6 7785 us/call
        db_displace       1.124 ms 674.6 kcyc   2.3% x6 187 us/call
  canvas_clear            90.98 us  54.6 kcyc   0.2% x1 91.0 us/call
  canvas_buffer_wait       1.03 us    621 cyc   0.0% x1 1.03 us/call
```

Wall min/mean/max was 45.884/48.502/50.929 ms. Render averaged 48.502 ms per frame, including 48.250 ms in `db_timeline_step`.

This regime leaves 13.998 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 2849–2864)

```
frame                    46.190 ms 27.71 Mcyc 100.0%
  pov_preserve_half      142.86 us  85.7 kcyc   0.3% x1 143 us/call
  db_timeline_step       45.937 ms 27.56 Mcyc  99.5% x1 45937 us/call
    db_draw              45.820 ms 27.49 Mcyc  99.2% x1 45820 us/call
      db_draw_scene      45.806 ms 27.48 Mcyc  99.2% x1 45806 us/call
        db_mesh_plot     44.428 ms 26.66 Mcyc  96.2% x6 7405 us/call
        db_displace       1.126 ms 675.4 kcyc   2.4% x6 188 us/call
  canvas_clear            91.24 us  54.7 kcyc   0.2% x1 91.2 us/call
  canvas_buffer_wait       1.03 us    621 cyc   0.0% x1 1.03 us/call
```

Wall min/mean/max was 41.685/46.189/48.415 ms. Render averaged 46.189 ms per frame, including 45.937 ms in `db_timeline_step`.

This regime leaves 16.311 ms of mean render headroom per display window.

### Final complete in-show window (frames 3169–3184)

```
frame                    30.766 ms 18.46 Mcyc 100.0%
  pov_preserve_half      143.92 us  86.4 kcyc   0.5% x1 144 us/call
  db_timeline_step       30.514 ms 18.31 Mcyc  99.2% x1 30514 us/call
    db_draw              30.399 ms 18.24 Mcyc  98.8% x1 30399 us/call
      db_draw_scene      30.385 ms 18.23 Mcyc  98.8% x1 30385 us/call
        db_mesh_plot     29.668 ms 17.80 Mcyc  96.4% x6 4945 us/call
        db_displace      577.84 us 346.7 kcyc   1.9% x6 96.3 us/call
  canvas_clear            89.86 us  53.9 kcyc   0.3% x1 89.9 us/call
  canvas_buffer_wait       0.76 us    458 cyc   0.0% x1 0.76 us/call
```

Wall min/mean/max was 29.825/30.766/32.493 ms. Render averaged 30.765 ms per frame, including 30.514 ms in `db_timeline_step`.

This regime leaves 31.735 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| 9 | held preset | 15 | 2769–2784 | — | — | 48.250 | 48.502 | 51.441 | 16 | 0/320 |
| 1 | held preset | 15 | 209–224 | — | — | 44.400 | 44.653 | 47.355 | 16 | 0/318 |
| 10 | held preset | 15 | 3089–3104 | — | — | 35.844 | 36.096 | 39.854 | 16 | 0/309 |
| 8 | held preset | 15 | 2417–2432 | — | — | 35.412 | 35.664 | 38.120 | 16 | 0/320 |
| 4 | held preset | 15 | 1137–1152 | — | — | 27.517 | 27.770 | 31.408 | 16 | 0/320 |
| 3 | held preset | 15 | 849–864 | — | — | 25.828 | 26.080 | 28.510 | 16 | 0/320 |
| 2 | held preset | 15 | 353–368 | — | — | 24.866 | 25.118 | 28.131 | 16 | 0/320 |
| 7 | held preset | 15 | 1953–1968 | — | — | 23.484 | 23.733 | 25.602 | 16 | 0/320 |
| 6 | held preset | 15 | 1633–1648 | — | — | 17.261 | 17.512 | 20.130 | 16 | 0/320 |
| 5 | held preset | 15 | 1473–1488 | — | — | 16.643 | 16.896 | 18.330 | 16 | 0/320 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 2737–2752, measured over the separate 1,006,121 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1159.1/frame  0.500/1.638/20.573 us  3.02%
  isr_pack           144.9/frame  5.115/6.416/8.893 us  1.47%
  isr_dma_submit     144.9/frame  0.706/0.918/4.413 us  0.21%
```

Inclusive wake time leaves 60.612 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `db_timeline_step` took 47.883 ms/frame (99.5% of draw time).
2. `db_draw` took 47.737 ms/frame (99.2% of draw time).
3. `db_draw_scene` took 47.722 ms/frame (99.1% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared Plot raster, column-walk and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=DreamBalls` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
