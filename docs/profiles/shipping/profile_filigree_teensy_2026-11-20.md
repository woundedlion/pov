# Filigree on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 105-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/filigree_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Filigree at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 105 s slot, 111 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Filigree profile 111 16` from the captured revision |

Captured 2026-11-20 15:45, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `4ad9b64980ffdc9de2d8b0f62a61f03b2e93955267409cdfcc2de4ef94a06b95`.

Profile image size `FLASH: code:116608, data:456380, headers:8644 free for files:1449984 / RAM1: variables:314144, code:28632, padding:4136 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `73e82e585d36a80aec89958d13bc561aa4994e3b058108e453674e905343c39c`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 769–784. The root counter measured 852682.972 µs, compared with 852,681 µs of draw time (2.31 ppm).

## Frame cadence

**Peak confirmed live render 55.155 ms** at frame 1080 (jade-weave → ember-lattice). **0/1,666 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 78.329 ms. Live render averaged 42.602 ms.

Among complete live readouts, the largest window mean was 53.291 ms/frame at frames 785–800. The dominant scope averaged 41.791 ms/frame across those readouts.

The shipping ELF assigns 840 revolutions (105 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 105 s. Frames 2–1667 are confirmed live. The raw file contains 1,752 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 13 frame(s), 1668–1680, with uncertain membership. Their largest render was 45.085 ms and 0 exceeded the budget. The live denominator lies between 1,666 and 1,679. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 55.155 (2), spilled 🟢 0/1666 (0.00%).

README setup frame 1: 78.329 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    46.889 ms 28.13 Mcyc 100.0%
  pov_preserve_half      125.21 us  75.1 kcyc   0.3% x0.94 134 us/call
  fx_shader_draw         45.580 ms 27.35 Mcyc  97.2% x1 45580 us/call
    fx_prepare_column     5.668 ms  3.40 Mcyc  12.1% x144 39.4 us/call
  fx_palette_base_bake   712.24 us 427.3 kcyc   1.5% x1 712 us/call
    fx_palette_rgb_bake  635.08 us 381.0 kcyc   1.4% x1 635 us/call
  fx_prepare_frame        51.32 us  30.8 kcyc   0.1% x1 51.3 us/call
  fx_advance             175.22 us 105.1 kcyc   0.4% x1 175 us/call
  fx_timeline_step         7.24 us   4.3 kcyc   0.0% x1 7.24 us/call
  canvas_clear           104.03 us  62.4 kcyc   0.2% x1 104 us/call
  canvas_buffer_wait       1.82 us   1.1 kcyc   0.0% x1 1.82 us/call
```

Wall min/mean/max was 44.483/46.889/78.332 ms. Render averaged 46.888 ms per frame, including 45.580 ms in `fx_shader_draw`.

This regime leaves 15.612 ms of mean render headroom per display window.

### Window containing the live peak (frames 1073–1088)

```
frame                    41.915 ms 25.15 Mcyc 100.0%
  pov_preserve_half      134.65 us  80.8 kcyc   0.3% x1 135 us/call
  fx_shader_draw         41.106 ms 24.66 Mcyc  98.1% x1 41106 us/call
    fx_prepare_column     3.862 ms  2.32 Mcyc   9.2% x144 26.8 us/call
  fx_palette_base_bake    70.07 us  42.0 kcyc   0.2% x0.94 74.7 us/call
  fx_prepare_frame       231.77 us 139.1 kcyc   0.6% x1 232 us/call
  fx_advance              50.03 us  30.0 kcyc   0.1% x1 50.0 us/call
  fx_timeline_step        63.97 us  38.4 kcyc   0.2% x1 64.0 us/call
  canvas_clear            95.08 us  57.0 kcyc   0.2% x1 95.1 us/call
  canvas_buffer_wait       3.72 us   2.2 kcyc   0.0% x1 3.72 us/call
```

Wall min/mean/max was 31.964/41.915/55.156 ms. Render averaged 41.911 ms per frame, including 41.106 ms in `fx_shader_draw`.

This regime leaves 20.589 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 81–96)

```
frame                    44.989 ms 26.99 Mcyc 100.0%
  pov_preserve_half      133.32 us  80.0 kcyc   0.3% x1 133 us/call
  fx_shader_draw         43.854 ms 26.31 Mcyc  97.5% x1 43854 us/call
    fx_prepare_column     5.648 ms  3.39 Mcyc  12.6% x144 39.2 us/call
  fx_palette_base_bake   672.87 us 403.7 kcyc   1.5% x1 673 us/call
    fx_palette_rgb_bake  599.54 us 359.7 kcyc   1.3% x1 600 us/call
  fx_prepare_frame        49.10 us  29.5 kcyc   0.1% x1 49.1 us/call
  fx_advance              51.74 us  31.0 kcyc   0.1% x1 51.7 us/call
  fx_timeline_step         4.88 us   2.9 kcyc   0.0% x1 4.88 us/call
  canvas_clear            95.14 us  57.1 kcyc   0.2% x1 95.1 us/call
  canvas_buffer_wait       1.28 us    767 cyc   0.0% x1 1.28 us/call
```

Wall min/mean/max was 44.552/44.988/45.829 ms. Render averaged 44.987 ms per frame, including 43.854 ms in `fx_shader_draw`.

This regime leaves 17.513 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 785–800)

```
frame                    53.293 ms 31.98 Mcyc 100.0%
  pov_preserve_half      134.45 us  80.7 kcyc   0.3% x1 134 us/call
  fx_shader_draw         52.603 ms 31.56 Mcyc  98.7% x1 52603 us/call
    fx_prepare_column     5.122 ms  3.07 Mcyc   9.6% x144 35.6 us/call
  fx_palette_base_bake    79.53 us  47.7 kcyc   0.1% x1 79.5 us/call
  fx_prepare_frame        49.57 us  29.7 kcyc   0.1% x1 49.6 us/call
  fx_advance              43.37 us  26.0 kcyc   0.1% x1 43.4 us/call
  fx_timeline_step       123.75 us  74.2 kcyc   0.2% x1 124 us/call
  canvas_clear            94.09 us  56.5 kcyc   0.2% x1 94.1 us/call
  canvas_buffer_wait       2.13 us   1.3 kcyc   0.0% x1 2.13 us/call
```

Wall min/mean/max was 51.884/53.292/55.123 ms. Render averaged 53.291 ms per frame, including 52.603 ms in `fx_shader_draw`.

This regime leaves 9.209 ms of mean render headroom per display window.

### Final complete in-show window (frames 1649–1664)

```
frame                    32.435 ms 19.46 Mcyc 100.0%
  pov_preserve_half      134.89 us  80.9 kcyc   0.4% x1 135 us/call
  fx_shader_draw         31.865 ms 19.12 Mcyc  98.2% x1 31865 us/call
    fx_prepare_column     1.836 ms  1.10 Mcyc   5.7% x144 12.8 us/call
  fx_palette_base_bake    69.50 us  41.7 kcyc   0.2% x1 69.5 us/call
  fx_prepare_frame        54.34 us  32.6 kcyc   0.2% x1 54.3 us/call
  fx_advance              51.64 us  31.0 kcyc   0.2% x1 51.6 us/call
  fx_timeline_step         4.87 us   2.9 kcyc   0.0% x1 4.87 us/call
  canvas_clear            93.15 us  55.9 kcyc   0.3% x1 93.2 us/call
  canvas_buffer_wait       3.70 us   2.2 kcyc   0.0% x1 3.70 us/call
```

Wall min/mean/max was 31.804/32.434/32.978 ms. Render averaged 32.431 ms per frame, including 31.865 ms in `fx_shader_draw`.

This regime leaves 30.069 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| jade-weave | held preset | 36 | 81–96 | — | — | 43.854 | 44.987 | 45.828 | 16 | 0/598 |
| ember-lattice | held preset | 36 | 1489–1504 | — | — | 32.253 | 32.814 | 55.155 | 16 | 0/1068 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 1073–1088, measured over the separate 981,899 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1131.2/frame  0.565/1.648/27.085 us  3.04%
  isr_pack           141.4/frame  4.910/6.283/8.948 us  1.44%
  isr_dma_submit     141.4/frame  0.626/0.918/5.315 us  0.21%
```

Inclusive wake time leaves 60.600 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 41.106 ms/frame (98.1% of draw time).
2. `fx_prepare_column` took 3.862 ms/frame (9.2% of draw time).
3. `fx_prepare_frame` took 0.232 ms/frame (0.6% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Filigree` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
