# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/cosmiceyeball_ship.log`. Supersedes `profile_cosmiceyeball_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | CosmicEyeball at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh CosmicEyeball profile 66 32` from the captured revision |

Captured 2026-11-20 15:55, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `04c35d49674d7b24fe077c7e4d8f5da824ad0435ead930f9b2bc884fb6af90be`.

Profile image size `FLASH: code:97312, data:152024, headers:8712 free for files:1773568 / RAM1: variables:314144, code:26536, padding:6232 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `1f314e6b9d12e640c9790161296291f88625cc5ddb93623acd4d877363dc038a`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 385–416. The root counter measured 815770.630 µs, compared with 815,773 µs of draw time (2.91 ppm).

## Frame cadence

**Peak confirmed live render 27.513 ms** at frame 371 (CosmicEyeball). **0/943 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 47.003 ms. Live render averaged 25.260 ms.

Among complete live readouts, the largest window mean was 25.591 ms/frame at frames 865–896. The dominant scope averaged 22.528 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–944 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 945–960, with uncertain membership. Their largest render was 26.951 ms and 0 exceeded the budget. The live denominator lies between 943 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 27.513, spilled 🟢 0/943 (0.00%).

README setup frame 1: 47.003 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    25.709 ms 15.43 Mcyc 100.0%
  pov_preserve_half      139.62 us  83.8 kcyc   0.5% x0.97 144 us/call
  fx_shader_draw         23.092 ms 13.86 Mcyc  89.8% x1 23092 us/call
    fx_prepare_column    599.40 us 359.6 kcyc   2.3% x144 4.16 us/call
  fx_palette_base_bake   566.59 us 340.0 kcyc   2.2% x1 567 us/call
    fx_palette_rgb_bake  437.43 us 262.5 kcyc   1.7% x1 437 us/call
  fx_prepare_frame       478.63 us 287.2 kcyc   1.9% x1 479 us/call
  fx_advance              1.055 ms 632.8 kcyc   4.1% x1 1055 us/call
    palette_cycle_step    1.015 ms 608.9 kcyc   3.9% x1 1015 us/call
  fx_timeline_step        99.85 us  59.9 kcyc   0.4% x1 99.8 us/call
  canvas_clear            99.43 us  59.7 kcyc   0.4% x1 99.4 us/call
  canvas_buffer_wait       1.46 us    877 cyc   0.0% x1 1.46 us/call
```

Wall min/mean/max was 23.518/25.709/47.003 ms. Render averaged 25.708 ms per frame, including 23.092 ms in `fx_shader_draw`.

This regime leaves 36.792 ms of mean render headroom per display window.

### Window containing the live peak (frames 353–384)

```
frame                    25.414 ms 15.25 Mcyc 100.0%
  pov_preserve_half      143.37 us  86.0 kcyc   0.6% x1 143 us/call
  fx_shader_draw         22.665 ms 13.60 Mcyc  89.2% x1 22665 us/call
    fx_prepare_column    603.39 us 362.0 kcyc   2.4% x144 4.19 us/call
  fx_palette_base_bake   744.95 us 447.0 kcyc   2.9% x1 745 us/call
    fx_palette_rgb_bake  609.61 us 365.8 kcyc   2.4% x1 610 us/call
  fx_prepare_frame       481.44 us 288.9 kcyc   1.9% x1 481 us/call
  fx_advance              1.009 ms 605.2 kcyc   4.0% x1 1009 us/call
    palette_cycle_step   971.26 us 582.8 kcyc   3.8% x1 971 us/call
  fx_timeline_step        97.23 us  58.3 kcyc   0.4% x1 97.2 us/call
  canvas_clear            94.69 us  56.8 kcyc   0.4% x1 94.7 us/call
  canvas_buffer_wait       1.80 us   1.1 kcyc   0.0% x1 1.80 us/call
```

Wall min/mean/max was 23.563/25.414/27.515 ms. Render averaged 25.412 ms per frame, including 22.665 ms in `fx_shader_draw`.

This regime leaves 37.088 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    25.430 ms 15.26 Mcyc 100.0%
  pov_preserve_half      147.13 us  88.3 kcyc   0.6% x1 147 us/call
  fx_shader_draw         22.462 ms 13.48 Mcyc  88.3% x1 22462 us/call
    fx_prepare_column    597.02 us 358.2 kcyc   2.3% x144 4.15 us/call
  fx_palette_base_bake   908.09 us 544.9 kcyc   3.6% x1 908 us/call
    fx_palette_rgb_bake  775.19 us 465.1 kcyc   3.0% x1 775 us/call
  fx_prepare_frame       479.05 us 287.4 kcyc   1.9% x1 479 us/call
  fx_advance              1.056 ms 633.4 kcyc   4.2% x1 1056 us/call
    palette_cycle_step    1.015 ms 608.7 kcyc   4.0% x1 1015 us/call
  fx_timeline_step        99.43 us  59.7 kcyc   0.4% x1 99.4 us/call
  canvas_clear            94.56 us  56.7 kcyc   0.4% x1 94.6 us/call
  canvas_buffer_wait       1.69 us   1.0 kcyc   0.0% x1 1.69 us/call
```

Wall min/mean/max was 23.811/25.430/27.142 ms. Render averaged 25.429 ms per frame, including 22.462 ms in `fx_shader_draw`.

This regime leaves 37.071 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 353–384, measured over the separate 2,000,510 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.3/frame  0.576/1.783/23.746 us  3.28%
  isr_pack           144.1/frame  6.986/7.738/10.263 us  1.78%
  isr_dma_submit     144.1/frame  0.596/0.918/9.228 us  0.21%
```

Inclusive wake time leaves 60.450 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 22.665 ms/frame (89.2% of draw time).
2. `fx_advance` took 1.009 ms/frame (4.0% of draw time).
3. `palette_cycle_step` took 0.971 ms/frame (3.8% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=CosmicEyeball` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
