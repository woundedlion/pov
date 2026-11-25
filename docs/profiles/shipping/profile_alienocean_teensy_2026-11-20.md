# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 30-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/alienocean_ship.log`. Supersedes `profile_alienocean_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | AlienOcean at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 30 s slot, 36 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh AlienOcean profile 36 32` from the captured revision |

Captured 2026-11-20 16:05, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `0e47265d9049293cc2ef6f0ff3305e0847c355d315ac514e3a2ef43f3fd9633c`.

Profile image size `FLASH: code:99976, data:152604, headers:8540 free for files:1770496 / RAM1: variables:314144, code:27800, padding:4968 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `4e0a70f5b29de3f4874f18ac739d32724a934f60cdd0bc2f9c89e989c1eef77c`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 321–352. The root counter measured 848482.973 µs, compared with 848,484 µs of draw time (1.21 ppm).

## Frame cadence

**Peak confirmed live render 28.276 ms** at frame 395 (AlienOcean). **0/463 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 51.217 ms. Live render averaged 26.411 ms.

Among complete live readouts, the largest window mean was 26.590 ms/frame at frames 385–416. The dominant scope averaged 23.783 ms/frame across those readouts.

The shipping ELF assigns 240 revolutions (30 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 30 s. Frames 2–464 are confirmed live. The raw file contains 553 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 465–480, with uncertain membership. Their largest render was 28.154 ms and 0 exceeded the budget. The live denominator lies between 463 and 479. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 28.276, spilled 🟢 0/463 (0.00%).

README setup frame 1: 51.217 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    26.961 ms 16.18 Mcyc 100.0%
  pov_preserve_half      140.00 us  84.0 kcyc   0.5% x0.97 145 us/call
  fx_shader_draw         24.415 ms 14.65 Mcyc  90.6% x1 24415 us/call
    fx_prepare_column    723.50 us 434.1 kcyc   2.7% x144 5.02 us/call
  fx_palette_base_bake   533.76 us 320.3 kcyc   2.0% x1 534 us/call
    fx_palette_rgb_bake  432.71 us 259.6 kcyc   1.6% x1 433 us/call
  fx_prepare_frame       599.07 us 359.4 kcyc   2.2% x1 599 us/call
  fx_advance             913.02 us 547.8 kcyc   3.4% x1 913 us/call
    palette_cycle_step   882.16 us 529.3 kcyc   3.3% x1 882 us/call
  fx_timeline_step        62.57 us  37.5 kcyc   0.2% x1 62.6 us/call
  canvas_clear            99.03 us  59.4 kcyc   0.4% x1 99.0 us/call
  canvas_buffer_wait       3.36 us   2.0 kcyc   0.0% x1 3.36 us/call
```

Wall min/mean/max was 24.671/26.960/51.217 ms. Render averaged 26.957 ms per frame, including 24.415 ms in `fx_shader_draw`.

This regime leaves 35.543 ms of mean render headroom per display window.

### Window containing the live peak (frames 385–416)

```
frame                    26.592 ms 15.95 Mcyc 100.0%
  pov_preserve_half      147.39 us  88.4 kcyc   0.6% x1 147 us/call
  fx_shader_draw         23.743 ms 14.25 Mcyc  89.3% x1 23743 us/call
    fx_prepare_column    722.80 us 433.7 kcyc   2.7% x144 5.02 us/call
  fx_palette_base_bake   870.88 us 522.5 kcyc   3.3% x1 871 us/call
    fx_palette_rgb_bake  763.00 us 457.8 kcyc   2.9% x1 763 us/call
  fx_prepare_frame       561.69 us 337.0 kcyc   2.1% x1 562 us/call
  fx_advance             920.21 us 552.1 kcyc   3.5% x1 920 us/call
    palette_cycle_step   887.91 us 532.7 kcyc   3.3% x1 888 us/call
  fx_timeline_step        60.66 us  36.4 kcyc   0.2% x1 60.7 us/call
  canvas_clear            95.17 us  57.1 kcyc   0.4% x1 95.2 us/call
  canvas_buffer_wait       1.58 us    948 cyc   0.0% x1 1.58 us/call
```

Wall min/mean/max was 25.159/26.591/28.277 ms. Render averaged 26.590 ms per frame, including 23.743 ms in `fx_shader_draw`.

This regime leaves 35.910 ms of mean render headroom per display window.

### Final complete in-show window (frames 417–448)

```
frame                    26.575 ms 15.95 Mcyc 100.0%
  pov_preserve_half      145.62 us  87.4 kcyc   0.5% x1 146 us/call
  fx_shader_draw         23.757 ms 14.25 Mcyc  89.4% x1 23757 us/call
    fx_prepare_column    731.91 us 439.1 kcyc   2.8% x144 5.08 us/call
  fx_palette_base_bake   823.51 us 494.1 kcyc   3.1% x1 824 us/call
    fx_palette_rgb_bake  714.21 us 428.5 kcyc   2.7% x1 714 us/call
  fx_prepare_frame       555.53 us 333.3 kcyc   2.1% x1 556 us/call
  fx_advance             941.65 us 565.0 kcyc   3.5% x1 942 us/call
    palette_cycle_step   907.16 us 544.3 kcyc   3.4% x1 907 us/call
  fx_timeline_step        62.15 us  37.3 kcyc   0.2% x1 62.2 us/call
  canvas_clear            94.98 us  57.0 kcyc   0.4% x1 95.0 us/call
  canvas_buffer_wait       2.18 us   1.3 kcyc   0.0% x1 2.18 us/call
```

Wall min/mean/max was 25.128/26.575/28.176 ms. Render averaged 26.573 ms per frame, including 23.757 ms in `fx_shader_draw`.

This regime leaves 35.927 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 385–416, measured over the separate 2,000,046 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.1/frame  0.493/1.755/14.580 us  3.23%
  isr_pack           144.0/frame  5.303/7.631/10.326 us  1.75%
  isr_dma_submit     144.0/frame  0.591/0.916/3.520 us  0.21%
```

Inclusive wake time leaves 60.481 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 23.743 ms/frame (89.3% of draw time).
2. `fx_advance` took 0.920 ms/frame (3.5% of draw time).
3. `palette_cycle_step` took 0.888 ms/frame (3.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=AlienOcean` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
