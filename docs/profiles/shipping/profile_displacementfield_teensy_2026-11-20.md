# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 121-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/displacementfield_ship.log`. Supersedes `profile_displacementfield_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | DisplacementField at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 121 s slot, 127 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh DisplacementField profile 127 16` from the captured revision |

Captured 2026-11-20 15:17, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `7646adefe24a043d92052884e8a61fe380eb74506fd35c7861b39f7d006f56bb`.

Profile image size `FLASH: code:118136, data:150736, headers:8632 free for files:1754112 / RAM1: variables:314464, code:25272, padding:7496 free for local variables:177056 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `788ecea2f86b002e4c70e5624c0450ed9f49bf8187e4c8e9a7a78d8ab2648e37`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 593–608. The root counter measured 695612.532 µs, compared with 695,613 µs of draw time (0.67 ppm).

## Frame cadence

**Peak confirmed live render 44.765 ms** at frame 603 (noise). **0/1,921 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 18.584 ms. Live render averaged 35.847 ms.

Among complete live readouts, the largest window mean was 43.473 ms/frame at frames 593–608. The dominant scope averaged 35.648 ms/frame across those readouts.

The shipping ELF assigns 968 revolutions (121 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 121 s. Frames 2–1922 are confirmed live. The raw file contains 2,008 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 1923–1936, with uncertain membership. Their largest render was 30.819 ms and 0 exceeded the budget. The live denominator lies between 1,921 and 1,935. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 44.765 (2), spilled 🟢 0/1921 (0.00%).

README setup frame 1: 18.584 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    26.269 ms 15.76 Mcyc 100.0%
  pov_preserve_half      137.51 us  82.5 kcyc   0.5% x0.94 147 us/call
  df_timeline_step       26.011 ms 15.61 Mcyc  99.0% x1 26011 us/call
    df_draw_rings        25.891 ms 15.53 Mcyc  98.6% x1 25891 us/call
      df_fused_scan      10.003 ms  6.00 Mcyc  38.1% x0.56 17782 us/call
      df_lut_bake         6.860 ms  4.12 Mcyc  26.1% x26.81 256 us/call
        df_hue_table_prep 194.97 us 117.0 kcyc   0.7% x8.06 24.2 us/call
        df_octave_noise   4.726 ms  2.84 Mcyc  18.0% x26.81 176 us/call
      df_chunk_cull      294.92 us 176.9 kcyc   1.1% x26.88 11.0 us/call
      df_lut_bake         1.510 ms 905.9 kcyc   5.7% x18 83.9 us/call
        df_hue_table_prep  25.79 us  15.5 kcyc   0.1% x5.25 4.91 us/call
      df_chunk_cull      347.59 us 208.6 kcyc   1.3% x18 19.3 us/call
      df_fused_scan       4.897 ms  2.94 Mcyc  18.6% x0.44 11192 us/call
        ring_stack_table  1.629 ms 977.2 kcyc   6.2% x1.5 1086 us/call
        ring_stack_table  87.96 us  52.8 kcyc   0.3% x0.06 1407 us/call
    animation_noise_clock   0.17 us    104 cyc   0.0% x1 0.17 us/call
  canvas_clear            98.46 us  59.1 kcyc   0.4% x1 98.5 us/call
  canvas_buffer_wait       2.92 us   1.7 kcyc   0.0% x1 2.92 us/call
  df_prepare_fields        1.23 us    737 cyc   0.0% x1 1.23 us/call
```

Wall min/mean/max was 18.001/26.268/32.713 ms. Render averaged 26.266 ms per frame, including 26.011 ms in `df_timeline_step`.

This regime leaves 36.234 ms of mean render headroom per display window.

### Window containing the live peak (frames 593–608)

```
ring_stack_table          1.902 ms  1.14 Mcyc   4.4% x2 951 us/call
frame                    43.476 ms 26.09 Mcyc 100.0%
  pov_preserve_half      144.43 us  86.7 kcyc   0.3% x1 144 us/call
  df_timeline_step       43.202 ms 25.92 Mcyc  99.4% x1 43202 us/call
    df_draw_rings        43.055 ms 25.83 Mcyc  99.0% x1 43055 us/call
      df_fused_scan      25.707 ms 15.42 Mcyc  59.1% x1 25707 us/call
      df_lut_bake        15.140 ms  9.08 Mcyc  34.8% x44.19 343 us/call
        df_hue_table_prep  1.132 ms 679.1 kcyc   2.6% x21.62 52.3 us/call
        df_octave_noise  10.425 ms  6.26 Mcyc  24.0% x44.19 236 us/call
      df_chunk_cull      508.38 us 305.0 kcyc   1.2% x45.88 11.1 us/call
    animation_noise_clock   0.17 us    102 cyc   0.0% x1 0.17 us/call
  canvas_clear            91.46 us  54.9 kcyc   0.2% x1 91.5 us/call
  canvas_buffer_wait       2.34 us   1.4 kcyc   0.0% x1 2.34 us/call
  df_prepare_fields        1.71 us   1.0 kcyc   0.0% x1 1.71 us/call
```

Wall min/mean/max was 41.722/43.475/44.767 ms. Render averaged 43.473 ms per frame, including 43.202 ms in `df_timeline_step`.

This regime leaves 19.027 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 737–752)

```
ring_stack_table          1.899 ms  1.14 Mcyc   4.4% x2 949 us/call
frame                    43.267 ms 25.96 Mcyc 100.0%
  pov_preserve_half      142.97 us  85.8 kcyc   0.3% x1 143 us/call
  df_timeline_step       42.994 ms 25.80 Mcyc  99.4% x1 42994 us/call
    df_draw_rings        42.855 ms 25.71 Mcyc  99.0% x1 42855 us/call
      df_fused_scan      25.762 ms 15.46 Mcyc  59.5% x1 25762 us/call
      df_lut_bake        14.913 ms  8.95 Mcyc  34.5% x43.5 343 us/call
        df_hue_table_prep  1.078 ms 646.7 kcyc   2.5% x22.19 48.6 us/call
        df_octave_noise  10.351 ms  6.21 Mcyc  23.9% x43.5 238 us/call
      df_chunk_cull      494.67 us 296.8 kcyc   1.1% x45.5 10.9 us/call
    animation_noise_clock   0.17 us    103 cyc   0.0% x1 0.17 us/call
  canvas_clear            92.80 us  55.7 kcyc   0.2% x1 92.8 us/call
  canvas_buffer_wait       2.41 us   1.4 kcyc   0.0% x1 2.41 us/call
  df_prepare_fields        1.54 us    924 cyc   0.0% x1 1.54 us/call
```

Wall min/mean/max was 42.372/43.267/44.172 ms. Render averaged 43.265 ms per frame, including 42.994 ms in `df_timeline_step`.

This regime leaves 19.235 ms of mean render headroom per display window.

### Final complete in-show window (frames 1905–1920)

```
frame                    15.732 ms  9.44 Mcyc 100.0%
  pov_preserve_half      146.57 us  87.9 kcyc   0.9% x1 147 us/call
  df_timeline_step       15.462 ms  9.28 Mcyc  98.3% x1 15462 us/call
    df_draw_rings        15.290 ms  9.17 Mcyc  97.2% x1 15290 us/call
      df_fused_scan       6.567 ms  3.94 Mcyc  41.7% x0.56 11674 us/call
      df_lut_bake         1.631 ms 978.7 kcyc  10.4% x5.56 293 us/call
        df_hue_table_prep  66.92 us  40.2 kcyc   0.4% x3.88 17.3 us/call
        df_octave_noise  901.42 us 540.8 kcyc   5.7% x4 225 us/call
      df_chunk_cull       74.17 us  44.5 kcyc   0.5% x6.5 11.4 us/call
      df_lut_bake         1.171 ms 702.8 kcyc   7.4% x12 97.6 us/call
        df_hue_table_prep  33.11 us  19.9 kcyc   0.2% x6.62 5.00 us/call
      df_chunk_cull      275.38 us 165.2 kcyc   1.8% x14.81 18.6 us/call
      df_fused_scan       4.233 ms  2.54 Mcyc  26.9% x0.44 9675 us/call
        ring_stack_table  1.180 ms 707.8 kcyc   7.5% x1.56 755 us/call
    animation_noise_clock   0.17 us    103 cyc   0.0% x1 0.17 us/call
  canvas_clear            94.69 us  56.8 kcyc   0.6% x1 94.7 us/call
  canvas_buffer_wait       2.24 us   1.3 kcyc   0.0% x1 2.24 us/call
  df_prepare_fields        5.22 us   3.1 kcyc   0.0% x1 5.22 us/call
```

Wall min/mean/max was 10.867/15.731/29.255 ms. Render averaged 15.729 ms per frame, including 15.462 ms in `df_timeline_step`.

This regime leaves 46.771 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| noise | full-gain noise | 36 | 593–608 | — | — | 43.202 | 43.473 | 44.765 | 16 | 0/909 |
| balls | ball spawning | 55 | 1617–1632 | — | — | 36.174 | 36.465 | 37.338 | 16 | 0/1012 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 593–608, measured over the separate 1,001,972 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1154.3/frame  0.503/1.728/17.465 us  3.18%
  isr_pack           144.3/frame  5.896/7.020/9.960 us  1.61%
  isr_dma_submit     144.3/frame  0.586/0.915/6.140 us  0.21%
```

Inclusive wake time leaves 60.512 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `df_timeline_step` took 43.202 ms/frame (99.4% of draw time).
2. `df_draw_rings` took 43.055 ms/frame (99.0% of draw time).
3. `df_fused_scan` took 25.707 ms/frame (59.1% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the effect-local ball/noise field evaluators and `core/render/scan/ring_stack.h`. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=DisplacementField` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
