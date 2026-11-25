# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 240-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/alienbrain_ship.log`. Supersedes `profile_alienbrain_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | AlienBrain at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 240 s slot, 246 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh AlienBrain profile 246 16` from the captured revision |

Captured 2026-11-20 15:05, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `7d87900b22a100a7c8d9f836f8a52660d81c1ff188a1586beb4c570659afdfb9`.

Profile image size `FLASH: code:102216, data:152852, headers:9124 free for files:1767424 / RAM1: variables:314144, code:27464, padding:5304 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `4f145b4270b94367ab115696361a627643389f5996cbab0c6a3f3abc63c1a845`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 1969–1984. The root counter measured 518320.235 µs, compared with 518,322 µs of draw time (3.41 ppm).

## Frame cadence

**Peak confirmed live render 37.554 ms** at frame 2042 (alien-brain-2 → alien-brain-3). **0/3,825 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 58.378 ms. Live render averaged 30.103 ms.

Among complete live readouts, the largest window mean was 33.035 ms/frame at frames 897–912. The dominant scope averaged 26.663 ms/frame across those readouts.

The shipping ELF assigns 1,920 revolutions (240 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 240 s. Frames 2–3826 are confirmed live. The raw file contains 3,912 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 3827–3840, with uncertain membership. Their largest render was 35.260 ms and 0 exceeded the budget. The live denominator lies between 3,825 and 3,839. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 37.554 (4), spilled 🟢 0/3825 (0.00%).

README setup frame 1: 58.378 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    31.235 ms 18.74 Mcyc 100.0%
  pov_preserve_half      131.32 us  78.8 kcyc   0.4% x0.94 140 us/call
  fx_shader_draw         28.173 ms 16.90 Mcyc  90.2% x1 28173 us/call
    fx_prepare_column    533.52 us 320.1 kcyc   1.7% x144 3.70 us/call
  fx_palette_base_bake   648.63 us 389.2 kcyc   2.1% x1 649 us/call
    fx_palette_rgb_bake  505.45 us 303.3 kcyc   1.6% x1 505 us/call
  fx_prepare_frame       747.87 us 448.7 kcyc   2.4% x1 748 us/call
  fx_advance              1.136 ms 681.5 kcyc   3.6% x1 1136 us/call
    palette_cycle_step    1.077 ms 646.4 kcyc   3.4% x1 1077 us/call
  fx_timeline_step        83.16 us  49.9 kcyc   0.3% x1 83.2 us/call
  canvas_clear           102.49 us  61.5 kcyc   0.3% x1 102 us/call
  canvas_buffer_wait       2.97 us   1.8 kcyc   0.0% x1 2.97 us/call
```

Wall min/mean/max was 27.580/31.234/58.383 ms. Render averaged 31.232 ms per frame, including 28.173 ms in `fx_shader_draw`.

This regime leaves 31.268 ms of mean render headroom per display window.

### Window containing the live peak (frames 2033–2048)

```
frame                    32.486 ms 19.49 Mcyc 100.0%
  pov_preserve_half      143.08 us  85.8 kcyc   0.4% x1 143 us/call
  fx_shader_draw         28.354 ms 17.01 Mcyc  87.3% x1 28354 us/call
    fx_prepare_column     1.081 ms 648.4 kcyc   3.3% x144 7.50 us/call
  fx_palette_base_bake    1.666 ms 999.5 kcyc   5.1% x1 1666 us/call
    fx_palette_rgb_bake   1.507 ms 904.3 kcyc   4.6% x1 1507 us/call
  fx_prepare_frame       746.66 us 448.0 kcyc   2.3% x1 747 us/call
  fx_advance              1.088 ms 652.6 kcyc   3.3% x1 1088 us/call
    palette_cycle_step    1.032 ms 619.2 kcyc   3.2% x1 1032 us/call
  fx_timeline_step       186.32 us 111.8 kcyc   0.6% x1 186 us/call
  canvas_clear            95.27 us  57.2 kcyc   0.3% x1 95.3 us/call
  canvas_buffer_wait       2.46 us   1.5 kcyc   0.0% x1 2.46 us/call
```

Wall min/mean/max was 29.936/32.486/37.556 ms. Render averaged 32.484 ms per frame, including 28.354 ms in `fx_shader_draw`.

This regime leaves 30.016 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 3521–3536)

```
frame                    31.335 ms 18.80 Mcyc 100.0%
  pov_preserve_half      140.91 us  84.5 kcyc   0.4% x1 141 us/call
  fx_shader_draw         27.944 ms 16.77 Mcyc  89.2% x1 27944 us/call
    fx_prepare_column    549.60 us 329.8 kcyc   1.8% x144 3.82 us/call
  fx_palette_base_bake    1.055 ms 633.2 kcyc   3.4% x1 1055 us/call
    fx_palette_rgb_bake  906.01 us 543.6 kcyc   2.9% x1 906 us/call
  fx_prepare_frame       732.71 us 439.6 kcyc   2.3% x1 733 us/call
  fx_advance              1.083 ms 650.0 kcyc   3.5% x1 1083 us/call
    palette_cycle_step    1.025 ms 615.2 kcyc   3.3% x1 1025 us/call
  fx_timeline_step        76.15 us  45.7 kcyc   0.2% x1 76.1 us/call
  canvas_clear            95.07 us  57.0 kcyc   0.3% x1 95.1 us/call
  canvas_buffer_wait       2.64 us   1.6 kcyc   0.0% x1 2.64 us/call
```

Wall min/mean/max was 29.464/31.334/33.267 ms. Render averaged 31.332 ms per frame, including 27.944 ms in `fx_shader_draw`.

This regime leaves 31.168 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 897–912)

```
frame                    33.037 ms 19.82 Mcyc 100.0%
  pov_preserve_half      142.63 us  85.6 kcyc   0.4% x1 143 us/call
  fx_shader_draw         28.517 ms 17.11 Mcyc  86.3% x1 28517 us/call
    fx_prepare_column     1.069 ms 641.3 kcyc   3.2% x144 7.42 us/call
  fx_palette_base_bake    2.130 ms  1.28 Mcyc   6.4% x1 2130 us/call
    fx_palette_rgb_bake   1.976 ms  1.19 Mcyc   6.0% x1 1976 us/call
  fx_prepare_frame       653.32 us 392.0 kcyc   2.0% x1 653 us/call
  fx_advance              1.103 ms 662.0 kcyc   3.3% x1 1103 us/call
    palette_cycle_step    1.048 ms 628.7 kcyc   3.2% x1 1048 us/call
  fx_timeline_step       180.99 us 108.6 kcyc   0.5% x1 181 us/call
  canvas_clear            96.08 us  57.7 kcyc   0.3% x1 96.1 us/call
  canvas_buffer_wait       2.73 us   1.6 kcyc   0.0% x1 2.73 us/call
```

Wall min/mean/max was 29.515/33.037/36.735 ms. Render averaged 33.035 ms per frame, including 28.517 ms in `fx_shader_draw`.

This regime leaves 29.465 ms of mean render headroom per display window.

### Final complete in-show window (frames 3809–3824)

```
frame                    30.720 ms 18.43 Mcyc 100.0%
  pov_preserve_half      139.03 us  83.4 kcyc   0.5% x1 139 us/call
  fx_shader_draw         26.522 ms 15.91 Mcyc  86.3% x1 26522 us/call
    fx_prepare_column    547.41 us 328.4 kcyc   1.8% x144 3.80 us/call
  fx_palette_base_bake    1.849 ms  1.11 Mcyc   6.0% x1 1849 us/call
    fx_palette_rgb_bake   1.694 ms  1.02 Mcyc   5.5% x1 1694 us/call
  fx_prepare_frame       704.98 us 423.0 kcyc   2.3% x1 705 us/call
  fx_advance              1.123 ms 674.0 kcyc   3.7% x1 1123 us/call
    palette_cycle_step    1.063 ms 637.7 kcyc   3.5% x1 1063 us/call
  fx_timeline_step        79.92 us  48.0 kcyc   0.3% x1 79.9 us/call
  canvas_clear            95.84 us  57.5 kcyc   0.3% x1 95.8 us/call
  canvas_buffer_wait       2.58 us   1.5 kcyc   0.0% x1 2.58 us/call
```

Wall min/mean/max was 27.903/30.719/35.429 ms. Render averaged 30.717 ms per frame, including 26.522 ms in `fx_shader_draw`.

This regime leaves 31.783 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| alien-brain-4 | held preset | 36 | 3521–3536 | — | — | 27.944 | 31.332 | 36.295 | 16 | 0/1069 |
| alien-brain | held preset | 36 | 321–336 | — | — | 27.754 | 31.366 | 36.637 | 16 | 0/598 |
| alien-brain-2 | held preset | 36 | 1441–1456 | — | — | 27.279 | 30.229 | 36.816 | 16 | 0/1079 |
| alien-brain-3 | held preset | 37 | 2593–2608 | — | — | 26.314 | 30.690 | 37.554 | 16 | 0/1079 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 2033–2048, measured over the separate 999,601 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.6/frame  0.500/1.698/18.118 us  3.13%
  isr_pack           143.9/frame  6.180/7.185/9.685 us  1.65%
  isr_dma_submit     143.9/frame  0.623/0.916/5.138 us  0.21%
```

Inclusive wake time leaves 60.544 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 28.354 ms/frame (87.3% of draw time).
2. `fx_palette_base_bake` took 1.666 ms/frame (5.1% of draw time).
3. `fx_palette_rgb_bake` took 1.507 ms/frame (4.6% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=AlienBrain` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
