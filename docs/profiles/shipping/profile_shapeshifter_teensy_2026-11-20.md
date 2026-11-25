# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 269-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/shapeshifter_ship.log`. Supersedes `profile_shapeshifter_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | ShapeShifter at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 269 s slot, 275 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh ShapeShifter profile 275 16` from the captured revision |

Captured 2026-11-20 15:33, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `d3a903e51c1fe5664235f29e06eaf9c21ad42a6d10c0e178322d6100fd8ed43a`.

Profile image size `FLASH: code:198104, data:152396, headers:8924 free for files:1672192 / RAM1: variables:314048, code:38728, padding:26808 free for local variables:144704 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `eef83e5595b9f0afd8a723ce7a7efd7a51fde3acd84c87e434e8add4dd618d19`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 3553–3568. The root counter measured 854257.137 µs, compared with 854,257 µs of draw time (0.16 ppm).

## Frame cadence

**Peak confirmed live render 56.029 ms** at frame 3554 (8). **0/4,289 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 76.265 ms. Live render averaged 34.517 ms.

Among complete live readouts, the largest window mean was 53.270 ms/frame at frames 3553–3568. The dominant scope averaged 34.132 ms/frame across those readouts.

The shipping ELF assigns 2,152 revolutions (269 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 269 s. Frames 2–4290 are confirmed live. The raw file contains 4,377 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 4291–4304, with uncertain membership. Their largest render was 44.014 ms and 0 exceeded the budget. The live denominator lies between 4,289 and 4,303. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 56.029 (9), spilled 🟢 0/4289 (0.00%).

README setup frame 1: 76.265 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    42.332 ms 25.40 Mcyc 100.0%
  pov_preserve_half      131.12 us  78.7 kcyc   0.3% x0.94 140 us/call
  ss_draw_all            41.936 ms 25.16 Mcyc  99.1% x1 41936 us/call
    ss_surface_knots      5.280 ms  3.17 Mcyc  12.5% x155.75 33.9 us/call
      ss_surface_mask     3.436 ms  2.06 Mcyc   8.1% x20.25 170 us/call
  ss_timeline_step        88.79 us  53.3 kcyc   0.2% x1 88.8 us/call
  ss_buffer_wait         128.73 us  77.2 kcyc   0.3% x1 129 us/call
    canvas_clear          99.41 us  59.6 kcyc   0.2% x1 99.4 us/call
    canvas_buffer_wait     1.35 us    809 cyc   0.0% x1 1.35 us/call
```

Wall min/mean/max was 38.783/42.332/76.267 ms. Render averaged 42.203 ms per frame, including 41.936 ms in `ss_draw_all`.

This regime leaves 20.297 ms of mean render headroom per display window.

### Window containing the live peak (frames 3553–3568)

```
frame                    53.391 ms 32.03 Mcyc 100.0%
  pov_preserve_half      134.89 us  80.9 kcyc   0.3% x1 135 us/call
  ss_draw_all            53.058 ms 31.83 Mcyc  99.4% x1 53058 us/call
  ss_timeline_step        65.09 us  39.1 kcyc   0.1% x1 65.1 us/call
  ss_buffer_wait         121.22 us  72.7 kcyc   0.2% x1 121 us/call
    canvas_clear          89.82 us  53.9 kcyc   0.2% x1 89.8 us/call
    canvas_buffer_wait     1.21 us    728 cyc   0.0% x1 1.21 us/call
```

Wall min/mean/max was 48.767/53.391/56.030 ms. Render averaged 53.270 ms per frame, including 53.058 ms in `ss_draw_all`.

This regime leaves 9.230 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 3329–3344)

```
frame                    45.730 ms 27.44 Mcyc 100.0%
  pov_preserve_half      136.56 us  81.9 kcyc   0.3% x1 137 us/call
  ss_draw_all            45.384 ms 27.23 Mcyc  99.2% x1 45384 us/call
  ss_timeline_step        71.89 us  43.1 kcyc   0.2% x1 71.9 us/call
  ss_buffer_wait         122.26 us  73.4 kcyc   0.3% x1 122 us/call
    canvas_clear          90.29 us  54.2 kcyc   0.2% x1 90.3 us/call
    canvas_buffer_wait     1.81 us   1.1 kcyc   0.0% x1 1.81 us/call
```

Wall min/mean/max was 44.149/45.729/47.910 ms. Render averaged 45.608 ms per frame, including 45.384 ms in `ss_draw_all`.

This regime leaves 16.892 ms of mean render headroom per display window.

### Final complete in-show window (frames 4273–4288)

```
frame                    28.630 ms 17.18 Mcyc 100.0%
  pov_preserve_half      142.44 us  85.5 kcyc   0.5% x1 142 us/call
  ss_draw_all            28.251 ms 16.95 Mcyc  98.7% x1 28251 us/call
  ss_timeline_step        81.00 us  48.6 kcyc   0.3% x1 81.0 us/call
  ss_buffer_wait         119.98 us  72.0 kcyc   0.4% x1 120 us/call
    canvas_clear          93.10 us  55.9 kcyc   0.3% x1 93.1 us/call
    canvas_buffer_wait     1.88 us   1.1 kcyc   0.0% x1 1.88 us/call
```

Wall min/mean/max was 27.704/28.629/29.597 ms. Render averaged 28.510 ms per frame, including 28.251 ms in `ss_draw_all`.

This regime leaves 33.990 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

Each solid has a 464-frame dwell and a 16-frame departure fade. The shape marker arrives at the fade midpoint. Qualified windows exclude both ends of the fade.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| 8 | held preset | 28 | 3553–3568 | — | — | 53.058 | 53.270 | 56.029 | 16 | 0/479 |
| 7 | held preset | 28 | 3137–3152 | — | — | 50.306 | 50.518 | 52.906 | 16 | 0/479 |
| 1 | held preset | 27 | 81–96 | — | — | 46.268 | 46.537 | 52.421 | 16 | 0/470 |
| 6 | held preset | 28 | 2721–2736 | — | — | 42.704 | 42.926 | 44.151 | 16 | 0/479 |
| 4 | held preset | 28 | 1809–1824 | — | — | 36.402 | 36.664 | 37.852 | 16 | 0/479 |
| 5 | held preset | 28 | 2145–2160 | — | — | 31.552 | 31.825 | 37.649 | 16 | 0/479 |
| 9 | held preset | 28 | 4097–4112 | — | — | 29.443 | 29.710 | 30.761 | 16 | 0/466 |
| 3 | held preset | 28 | 1217–1232 | — | — | 25.139 | 25.422 | 30.167 | 16 | 0/479 |
| 2 | held preset | 28 | 625–640 | — | — | 24.991 | 25.203 | 29.341 | 16 | 0/479 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 3553–3568, measured over the separate 998,056 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1149.8/frame  0.573/1.716/16.420 us  3.16%
  isr_pack           143.8/frame  6.300/7.140/9.620 us  1.64%
  isr_dma_submit     143.8/frame  0.625/0.920/8.953 us  0.21%
```

Inclusive wake time leaves 60.525 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `ss_draw_all` took 53.058 ms/frame (99.4% of draw time).
2. `pov_preserve_half` took 0.135 ms/frame (0.3% of draw time).
3. `canvas_clear` took 0.090 ms/frame (0.2% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared Plot raster, column-walk and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=ShapeShifter` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
