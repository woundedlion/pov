# KaleidoscopeHexBright on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 105-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/kaleidoscopehexbright_ship.log`. Supersedes `profile_kaleidoscopehexbright_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | KaleidoscopeHexBright at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 105 s slot, 111 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeHexBright profile 111 16` from the captured revision |

Captured 2026-11-20 15:22, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `3b128fc7346c338420a50b2c0a75477abacc25ceec29237d4f8ef74584a9556f`.

Profile image size `FLASH: code:100872, data:152304, headers:8968 free for files:1769472 / RAM1: variables:314144, code:26664, padding:6104 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `aa639ad4f7d18f2a3778a7b8fa1f91a0b02b08080bd62083dee86150716575f4`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 913–928. The root counter measured 594645.393 µs, compared with 594,645 µs of draw time (0.66 ppm).

## Frame cadence

**Peak confirmed live render 42.290 ms** at frame 395 (hex-twin-wave). **0/1,667 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 59.706 ms. Live render averaged 34.544 ms.

Among complete live readouts, the largest window mean was 37.749 ms/frame at frames 337–352. The dominant scope averaged 28.697 ms/frame across those readouts.

The shipping ELF assigns 840 revolutions (105 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 105 s. Frames 2–1668 are confirmed live. The raw file contains 1,752 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 12 frame(s), 1669–1680, with uncertain membership. Their largest render was 36.677 ms and 0 exceeded the budget. The live denominator lies between 1,667 and 1,679. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 42.290 (2), spilled 🟢 0/1667 (0.00%).

README setup frame 1: 59.706 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    33.958 ms 20.38 Mcyc 100.0%
  pov_preserve_half      133.11 us  79.9 kcyc   0.4% x0.94 142 us/call
  fx_shader_draw         28.278 ms 16.97 Mcyc  83.3% x1 28278 us/call
    fx_prepare_column    727.86 us 436.7 kcyc   2.1% x144 5.05 us/call
  fx_palette_base_bake   778.11 us 466.9 kcyc   2.3% x1 778 us/call
    fx_palette_rgb_bake  665.84 us 399.5 kcyc   2.0% x1 666 us/call
  fx_prepare_frame        3.340 ms  2.00 Mcyc   9.8% x1 3340 us/call
  fx_advance              1.076 ms 645.3 kcyc   3.2% x1 1076 us/call
    palette_cycle_step    1.022 ms 613.2 kcyc   3.0% x1 1022 us/call
  fx_timeline_step        86.83 us  52.1 kcyc   0.3% x1 86.8 us/call
  canvas_clear           100.99 us  60.6 kcyc   0.3% x1 101 us/call
  canvas_buffer_wait       1.61 us    967 cyc   0.0% x1 1.61 us/call
```

Wall min/mean/max was 29.911/33.958/59.708 ms. Render averaged 33.957 ms per frame, including 28.278 ms in `fx_shader_draw`.

This regime leaves 28.543 ms of mean render headroom per display window.

### Window containing the live peak (frames 385–400)

```
frame                    37.672 ms 22.60 Mcyc 100.0%
  pov_preserve_half      142.33 us  85.4 kcyc   0.4% x1 142 us/call
  fx_shader_draw         30.760 ms 18.46 Mcyc  81.7% x1 30760 us/call
    fx_prepare_column    765.05 us 459.0 kcyc   2.0% x144 5.31 us/call
  fx_palette_base_bake    2.282 ms  1.37 Mcyc   6.1% x1 2282 us/call
    fx_palette_rgb_bake   2.167 ms  1.30 Mcyc   5.8% x1 2167 us/call
  fx_prepare_frame        3.209 ms  1.93 Mcyc   8.5% x1 3209 us/call
  fx_advance             945.56 us 567.3 kcyc   2.5% x1 946 us/call
    palette_cycle_step   890.37 us 534.2 kcyc   2.4% x1 890 us/call
  fx_timeline_step        90.78 us  54.5 kcyc   0.2% x1 90.8 us/call
  canvas_clear            94.29 us  56.6 kcyc   0.3% x1 94.3 us/call
  canvas_buffer_wait       2.34 us   1.4 kcyc   0.0% x1 2.34 us/call
```

Wall min/mean/max was 33.334/37.671/42.293 ms. Render averaged 37.669 ms per frame, including 30.760 ms in `fx_shader_draw`.

This regime leaves 24.831 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 401–416)

```
frame                    37.594 ms 22.56 Mcyc 100.0%
  pov_preserve_half      140.02 us  84.0 kcyc   0.4% x1 140 us/call
  fx_shader_draw         31.015 ms 18.61 Mcyc  82.5% x1 31015 us/call
    fx_prepare_column    768.79 us 461.3 kcyc   2.0% x144 5.34 us/call
  fx_palette_base_bake    1.946 ms  1.17 Mcyc   5.2% x1 1946 us/call
    fx_palette_rgb_bake   1.829 ms  1.10 Mcyc   4.9% x1 1829 us/call
  fx_prepare_frame        3.178 ms  1.91 Mcyc   8.5% x1 3178 us/call
  fx_advance             970.60 us 582.4 kcyc   2.6% x1 971 us/call
    palette_cycle_step   917.11 us 550.3 kcyc   2.4% x1 917 us/call
  fx_timeline_step        93.67 us  56.2 kcyc   0.2% x1 93.7 us/call
  canvas_clear            93.66 us  56.2 kcyc   0.2% x1 93.7 us/call
  canvas_buffer_wait       1.76 us   1.1 kcyc   0.0% x1 1.76 us/call
```

Wall min/mean/max was 34.148/37.593/42.019 ms. Render averaged 37.592 ms per frame, including 31.015 ms in `fx_shader_draw`.

This regime leaves 24.908 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 913–928)

```
frame                    37.165 ms 22.30 Mcyc 100.0%
  pov_preserve_half      139.99 us  84.0 kcyc   0.4% x1 140 us/call
  fx_shader_draw         31.098 ms 18.66 Mcyc  83.7% x1 31098 us/call
    fx_prepare_column     1.133 ms 679.8 kcyc   3.0% x144 7.87 us/call
  fx_palette_base_bake    1.318 ms 790.7 kcyc   3.5% x1 1318 us/call
    fx_palette_rgb_bake   1.187 ms 712.0 kcyc   3.2% x1 1187 us/call
  fx_prepare_frame        3.251 ms  1.95 Mcyc   8.7% x1 3251 us/call
  fx_advance             918.40 us 551.0 kcyc   2.5% x1 918 us/call
    palette_cycle_step   864.99 us 519.0 kcyc   2.3% x1 865 us/call
  fx_timeline_step       180.46 us 108.3 kcyc   0.5% x1 180 us/call
  canvas_clear            93.67 us  56.2 kcyc   0.3% x1 93.7 us/call
  canvas_buffer_wait       1.88 us   1.1 kcyc   0.0% x1 1.88 us/call
```

Wall min/mean/max was 34.536/37.165/39.859 ms. Render averaged 37.163 ms per frame, including 31.098 ms in `fx_shader_draw`.

This regime leaves 25.337 ms of mean render headroom per display window.

### Final complete in-show window (frames 1649–1664)

```
frame                    35.524 ms 21.31 Mcyc 100.0%
  pov_preserve_half      140.27 us  84.2 kcyc   0.4% x1 140 us/call
  fx_shader_draw         29.293 ms 17.58 Mcyc  82.5% x1 29293 us/call
    fx_prepare_column    742.56 us 445.5 kcyc   2.1% x144 5.16 us/call
  fx_palette_base_bake    1.506 ms 903.6 kcyc   4.2% x1 1506 us/call
    fx_palette_rgb_bake   1.387 ms 832.2 kcyc   3.9% x1 1387 us/call
  fx_prepare_frame        3.374 ms  2.02 Mcyc   9.5% x1 3374 us/call
  fx_advance             875.76 us 525.5 kcyc   2.5% x1 876 us/call
    palette_cycle_step   824.69 us 494.8 kcyc   2.3% x1 825 us/call
  fx_timeline_step        87.47 us  52.5 kcyc   0.2% x1 87.5 us/call
  canvas_clear            94.12 us  56.5 kcyc   0.3% x1 94.1 us/call
  canvas_buffer_wait       1.92 us   1.1 kcyc   0.0% x1 1.92 us/call
```

Wall min/mean/max was 34.388/35.524/37.315 ms. Render averaged 35.522 ms per frame, including 29.293 ms in `fx_shader_draw`.

This regime leaves 26.978 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| hex-twin-wave | held preset | 36 | 401–416 | — | — | 31.015 | 37.592 | 42.290 | 16 | 0/598 |
| hex-twin-wave-alt | held preset | 36 | 1569–1584 | — | — | 29.943 | 36.099 | 39.922 | 16 | 0/1069 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 385–400, measured over the separate 1,002,131 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1154.6/frame  0.515/1.745/18.525 us  3.21%
  isr_pack           144.3/frame  6.735/7.590/10.003 us  1.74%
  isr_dma_submit     144.3/frame  0.680/0.916/6.001 us  0.21%
```

Inclusive wake time leaves 60.494 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 30.760 ms/frame (81.7% of draw time).
2. `fx_prepare_frame` took 3.209 ms/frame (8.5% of draw time).
3. `fx_palette_base_bake` took 2.282 ms/frame (6.1% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeHexBright` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
