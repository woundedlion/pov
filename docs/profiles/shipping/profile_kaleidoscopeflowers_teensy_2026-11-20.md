# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 173-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/kaleidoscopeflowers_ship.log`. Supersedes `profile_kaleidoscopeflowers_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | KaleidoscopeFlowers at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 173 s slot, 179 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeFlowers profile 179 16` from the captured revision |

Captured 2026-11-20 15:14, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `88c053b188d43b538306a864c6a20829bab490bdcbef841b466475417c72691f`.

Profile image size `FLASH: code:102672, data:153008, headers:8512 free for files:1767424 / RAM1: variables:314144, code:27784, padding:4984 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `d731a8ada035a1bbae933b6b18eabea24becc07290644c4e339941ed36a4be96`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 1953–1968. The root counter measured 648557.478 µs, compared with 648,559 µs of draw time (2.35 ppm).

## Frame cadence

**Peak confirmed live render 44.353 ms** at frame 1960 (open-grid → fine-grid). **0/2,755 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 63.854 ms. Live render averaged 34.803 ms.

Among complete live readouts, the largest window mean was 40.532 ms/frame at frames 1953–1968. The dominant scope averaged 31.433 ms/frame across those readouts.

The shipping ELF assigns 1,384 revolutions (173 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 173 s. Frames 2–2756 are confirmed live. The raw file contains 2,840 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 12 frame(s), 2757–2768, with uncertain membership. Their largest render was 36.952 ms and 0 exceeded the budget. The live denominator lies between 2,755 and 2,767. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 44.353 (3), spilled 🟢 0/2755 (0.00%).

README setup frame 1: 63.854 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    34.573 ms 20.74 Mcyc 100.0%
  pov_preserve_half      134.14 us  80.5 kcyc   0.4% x0.94 143 us/call
  fx_shader_draw         31.367 ms 18.82 Mcyc  90.7% x1 31367 us/call
    fx_prepare_column    841.57 us 504.9 kcyc   2.4% x144 5.84 us/call
  fx_palette_base_bake   796.72 us 478.0 kcyc   2.3% x1 797 us/call
    fx_palette_rgb_bake  661.21 us 396.7 kcyc   1.9% x1 661 us/call
  fx_prepare_frame       890.56 us 534.3 kcyc   2.6% x1 891 us/call
  fx_advance             986.52 us 591.9 kcyc   2.9% x1 987 us/call
    palette_cycle_step   936.70 us 562.0 kcyc   2.7% x1 937 us/call
  fx_timeline_step        85.52 us  51.3 kcyc   0.2% x1 85.5 us/call
  canvas_clear           104.07 us  62.4 kcyc   0.3% x1 104 us/call
  canvas_buffer_wait       2.15 us   1.3 kcyc   0.0% x1 2.15 us/call
```

Wall min/mean/max was 30.945/34.572/63.855 ms. Render averaged 34.570 ms per frame, including 31.367 ms in `fx_shader_draw`.

This regime leaves 27.930 ms of mean render headroom per display window.

### Window containing the live peak (frames 1953–1968)

```
frame                    40.535 ms 24.32 Mcyc 100.0%
  pov_preserve_half      141.94 us  85.2 kcyc   0.4% x1 142 us/call
  fx_shader_draw         35.805 ms 21.48 Mcyc  88.3% x1 35805 us/call
    fx_prepare_column     1.925 ms  1.16 Mcyc   4.7% x144 13.4 us/call
  fx_palette_base_bake    2.324 ms  1.39 Mcyc   5.7% x1 2324 us/call
    fx_palette_rgb_bake   2.179 ms  1.31 Mcyc   5.4% x1 2179 us/call
  fx_prepare_frame       816.84 us 490.1 kcyc   2.0% x1 817 us/call
  fx_advance             952.35 us 571.4 kcyc   2.3% x1 952 us/call
    palette_cycle_step   893.80 us 536.3 kcyc   2.2% x1 894 us/call
  fx_timeline_step       192.44 us 115.5 kcyc   0.5% x1 192 us/call
  canvas_clear            97.01 us  58.2 kcyc   0.2% x1 97.0 us/call
  canvas_buffer_wait       2.77 us   1.7 kcyc   0.0% x1 2.77 us/call
```

Wall min/mean/max was 37.861/40.534/44.356 ms. Render averaged 40.532 ms per frame, including 35.805 ms in `fx_shader_draw`.

This regime leaves 21.968 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 2641–2656)

```
frame                    37.191 ms 22.31 Mcyc 100.0%
  pov_preserve_half      144.02 us  86.4 kcyc   0.4% x1 144 us/call
  fx_shader_draw         33.754 ms 20.25 Mcyc  90.8% x1 33754 us/call
    fx_prepare_column    846.43 us 507.9 kcyc   2.3% x144 5.88 us/call
  fx_palette_base_bake    1.278 ms 766.7 kcyc   3.4% x1 1278 us/call
    fx_palette_rgb_bake   1.139 ms 683.3 kcyc   3.1% x1 1139 us/call
  fx_prepare_frame       700.15 us 420.1 kcyc   1.9% x1 700 us/call
  fx_advance             935.02 us 561.0 kcyc   2.5% x1 935 us/call
    palette_cycle_step   887.17 us 532.3 kcyc   2.4% x1 887 us/call
  fx_timeline_step        90.42 us  54.3 kcyc   0.2% x1 90.4 us/call
  canvas_clear            94.39 us  56.6 kcyc   0.3% x1 94.4 us/call
  canvas_buffer_wait       2.77 us   1.7 kcyc   0.0% x1 2.77 us/call
```

Wall min/mean/max was 35.461/37.190/38.671 ms. Render averaged 37.188 ms per frame, including 33.754 ms in `fx_shader_draw`.

This regime leaves 25.312 ms of mean render headroom per display window.

### Final complete in-show window (frames 2737–2752)

```
frame                    35.044 ms 21.03 Mcyc 100.0%
  pov_preserve_half      143.08 us  85.8 kcyc   0.4% x1 143 us/call
  fx_shader_draw         31.813 ms 19.09 Mcyc  90.8% x1 31813 us/call
    fx_prepare_column    845.99 us 507.6 kcyc   2.4% x144 5.87 us/call
  fx_palette_base_bake    1.129 ms 677.6 kcyc   3.2% x1 1129 us/call
    fx_palette_rgb_bake  988.44 us 593.1 kcyc   2.8% x1 988 us/call
  fx_prepare_frame       623.46 us 374.1 kcyc   1.8% x1 623 us/call
  fx_advance             947.90 us 568.7 kcyc   2.7% x1 948 us/call
    palette_cycle_step   897.83 us 538.7 kcyc   2.6% x1 898 us/call
  fx_timeline_step        93.32 us  56.0 kcyc   0.3% x1 93.3 us/call
  canvas_clear            94.97 us  57.0 kcyc   0.3% x1 95.0 us/call
  canvas_buffer_wait       2.97 us   1.8 kcyc   0.0% x1 2.97 us/call
```

Wall min/mean/max was 33.558/35.043/36.504 ms. Render averaged 35.041 ms per frame, including 31.813 ms in `fx_shader_draw`.

This regime leaves 27.459 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| fine-grid | held preset | 37 | 2641–2656 | — | — | 33.754 | 37.188 | 44.353 | 16 | 0/1078 |
| open-grid | held preset | 36 | 1633–1648 | — | — | 33.738 | 37.478 | 40.270 | 16 | 0/1079 |
| double-map | held preset | 36 | 193–208 | — | — | 32.600 | 36.376 | 40.180 | 16 | 0/598 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 1953–1968, measured over the separate 1,003,911 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1156.6/frame  0.583/1.766/25.773 us  3.25%
  isr_pack           144.6/frame  6.913/7.631/10.561 us  1.75%
  isr_dma_submit     144.6/frame  0.613/0.916/2.861 us  0.21%
```

Inclusive wake time leaves 60.469 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 35.805 ms/frame (88.3% of draw time).
2. `fx_palette_base_bake` took 2.324 ms/frame (5.7% of draw time).
3. `fx_palette_rgb_bake` took 2.179 ms/frame (5.4% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeFlowers` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
