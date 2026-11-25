# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 105-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/mobiusgrid_ship.log`. Supersedes `profile_mobiusgrid_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | MobiusGrid at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 105 s slot, 111 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh MobiusGrid profile 111 16` from the captured revision |

Captured 2026-11-20 15:49, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `77f191fbae7b4203e65be6ccf9dccf9f3994cd45f026723f4abe99a28073ccf3`.

Profile image size `FLASH: code:99304, data:152356, headers:8436 free for files:1771520 / RAM1: variables:314112, code:25016, padding:7752 free for local variables:177408 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `3f86cf4b9e6493bd4c7477f74210b8c868ca9ae6447cd457f415ef5e79a4d5d3`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 689–704. The root counter measured 522175.415 µs, compared with 522,174 µs of draw time (2.71 ppm).

## Frame cadence

**Peak confirmed live render 34.612 ms** at frame 700 (mobius-grid → mobius-grid-2). **0/1,666 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 59.131 ms. Live render averaged 31.497 ms.

Among complete live readouts, the largest window mean was 32.634 ms/frame at frames 689–704. The dominant scope averaged 29.128 ms/frame across those readouts.

The shipping ELF assigns 840 revolutions (105 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 105 s. Frames 2–1667 are confirmed live. The raw file contains 1,752 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 13 frame(s), 1668–1680, with uncertain membership. Their largest render was 32.997 ms and 0 exceeded the budget. The live denominator lies between 1,666 and 1,679. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 34.612 (2), spilled 🟢 0/1666 (0.00%).

README setup frame 1: 59.131 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    32.811 ms 19.69 Mcyc 100.0%
  pov_preserve_half      135.31 us  81.2 kcyc   0.4% x0.94 144 us/call
  fx_shader_draw         30.457 ms 18.27 Mcyc  92.8% x1 30457 us/call
    fx_prepare_column    569.86 us 341.9 kcyc   1.7% x144 3.96 us/call
  fx_palette_base_bake   137.50 us  82.5 kcyc   0.4% x1 138 us/call
  fx_prepare_frame       464.02 us 278.4 kcyc   1.4% x1 464 us/call
  fx_advance              1.218 ms 731.0 kcyc   3.7% x1 1218 us/call
    palette_cycle_step    1.150 ms 690.3 kcyc   3.5% x1 1150 us/call
  fx_timeline_step        86.94 us  52.2 kcyc   0.3% x1 86.9 us/call
  canvas_clear            97.28 us  58.4 kcyc   0.3% x1 97.3 us/call
  canvas_buffer_wait       2.36 us   1.4 kcyc   0.0% x1 2.36 us/call
```

Wall min/mean/max was 29.516/32.811/59.131 ms. Render averaged 32.809 ms per frame, including 30.457 ms in `fx_shader_draw`.

This regime leaves 29.691 ms of mean render headroom per display window.

### Window containing the live peak (frames 689–704)

```
frame                    32.636 ms 19.58 Mcyc 100.0%
  pov_preserve_half      141.64 us  85.0 kcyc   0.4% x1 142 us/call
  fx_shader_draw         30.153 ms 18.09 Mcyc  92.4% x1 30153 us/call
    fx_prepare_column     1.051 ms 630.6 kcyc   3.2% x144 7.30 us/call
  fx_palette_base_bake   142.39 us  85.4 kcyc   0.4% x1 142 us/call
  fx_prepare_frame       494.31 us 296.6 kcyc   1.5% x1 494 us/call
  fx_advance              1.211 ms 726.7 kcyc   3.7% x1 1211 us/call
    palette_cycle_step    1.139 ms 683.7 kcyc   3.5% x1 1139 us/call
  fx_timeline_step       181.32 us 108.8 kcyc   0.6% x1 181 us/call
  canvas_clear            88.94 us  53.4 kcyc   0.3% x1 88.9 us/call
  canvas_buffer_wait       2.33 us   1.4 kcyc   0.0% x1 2.33 us/call
```

Wall min/mean/max was 30.921/32.635/34.613 ms. Render averaged 32.634 ms per frame, including 30.153 ms in `fx_shader_draw`.

This regime leaves 29.866 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 577–592)

```
frame                    31.859 ms 19.12 Mcyc 100.0%
  pov_preserve_half      142.63 us  85.6 kcyc   0.4% x1 143 us/call
  fx_shader_draw         29.493 ms 17.70 Mcyc  92.6% x1 29493 us/call
    fx_prepare_column    622.41 us 373.4 kcyc   2.0% x144 4.32 us/call
  fx_palette_base_bake   151.67 us  91.0 kcyc   0.5% x1 152 us/call
  fx_prepare_frame       497.80 us 298.7 kcyc   1.6% x1 498 us/call
  fx_advance              1.183 ms 709.9 kcyc   3.7% x1 1183 us/call
    palette_cycle_step    1.112 ms 667.1 kcyc   3.5% x1 1112 us/call
  fx_timeline_step        94.32 us  56.6 kcyc   0.3% x1 94.3 us/call
  canvas_clear            88.97 us  53.4 kcyc   0.3% x1 89.0 us/call
  canvas_buffer_wait       2.41 us   1.4 kcyc   0.0% x1 2.41 us/call
```

Wall min/mean/max was 29.860/31.859/33.622 ms. Render averaged 31.857 ms per frame, including 29.493 ms in `fx_shader_draw`.

This regime leaves 30.643 ms of mean render headroom per display window.

### Final complete in-show window (frames 1649–1664)

```
frame                    31.190 ms 18.71 Mcyc 100.0%
  pov_preserve_half      141.86 us  85.1 kcyc   0.5% x1 142 us/call
  fx_shader_draw         28.856 ms 17.31 Mcyc  92.5% x1 28856 us/call
    fx_prepare_column    612.48 us 367.5 kcyc   2.0% x144 4.25 us/call
  fx_palette_base_bake   143.75 us  86.3 kcyc   0.5% x1 144 us/call
  fx_prepare_frame       485.75 us 291.4 kcyc   1.6% x1 486 us/call
  fx_advance              1.176 ms 705.5 kcyc   3.8% x1 1176 us/call
    palette_cycle_step    1.114 ms 668.4 kcyc   3.6% x1 1114 us/call
  fx_timeline_step        94.99 us  57.0 kcyc   0.3% x1 95.0 us/call
  canvas_clear            88.99 us  53.4 kcyc   0.3% x1 89.0 us/call
  canvas_buffer_wait       2.03 us   1.2 kcyc   0.0% x1 2.03 us/call
```

Wall min/mean/max was 29.483/31.190/32.906 ms. Render averaged 31.188 ms per frame, including 28.856 ms in `fx_shader_draw`.

This regime leaves 31.312 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| mobius-grid | held preset | 36 | 577–592 | — | — | 29.493 | 31.857 | 33.669 | 16 | 0/598 |
| mobius-grid-2 | held preset | 36 | 1377–1392 | — | — | 28.936 | 31.261 | 34.612 | 16 | 0/1068 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 689–704, measured over the separate 1,000,117 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.2/frame  0.495/1.735/16.591 us  3.20%
  isr_pack           144.1/frame  6.373/7.323/9.931 us  1.68%
  isr_dma_submit     144.1/frame  0.806/0.913/0.980 us  0.21%
```

Inclusive wake time leaves 60.500 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 30.153 ms/frame (92.4% of draw time).
2. `fx_advance` took 1.211 ms/frame (3.7% of draw time).
3. `palette_cycle_step` took 1.139 ms/frame (3.5% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=MobiusGrid` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
