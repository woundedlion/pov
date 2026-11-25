# KaleidoscopeMelt on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/kaleidoscopemelt_ship.log`. Supersedes `profile_kaleidoscopemelt_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | KaleidoscopeMelt at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh KaleidoscopeMelt profile 66 32` from the captured revision |

Captured 2026-11-20 16:18, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `3be93fd9b661578283e8dd54808bfb97ed1f3bb040f70ac588d4296ce746f186`.

Profile image size `FLASH: code:98336, data:151968, headers:8768 free for files:1772544 / RAM1: variables:314112, code:29432, padding:3336 free for local variables:177408 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `8584254c02dcccf459a341e0782dc3756a151512a5ce1318ed06e99419049e7a`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 193–224. The root counter measured 1575287.787 µs, compared with 1,575,285 µs of draw time (1.77 ppm).

## Frame cadence

**Peak confirmed live render 51.489 ms** at frame 191 (KaleidoscopeMelt). **0/950 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 81.763 ms. Live render averaged 47.367 ms.

Among complete live readouts, the largest window mean was 49.227 ms/frame at frames 193–224. The dominant scope averaged 43.100 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–951 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 9 frame(s), 952–960, with uncertain membership. Their largest render was 49.929 ms and 0 exceeded the budget. The live denominator lies between 950 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 51.489, spilled 🟢 0/950 (0.00%).

README setup frame 1: 81.763 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    47.108 ms 28.26 Mcyc 100.0%
  pov_preserve_half      134.80 us  80.9 kcyc   0.3% x0.97 139 us/call
  fx_shader_draw         42.799 ms 25.68 Mcyc  90.9% x1 42799 us/call
    fx_prepare_column    468.07 us 280.8 kcyc   1.0% x144 3.25 us/call
  fx_palette_base_bake   310.59 us 186.4 kcyc   0.7% x1 311 us/call
  fx_prepare_frame        2.437 ms  1.46 Mcyc   5.2% x1 2437 us/call
  fx_advance              1.074 ms 644.6 kcyc   2.3% x1 1074 us/call
    palette_cycle_step    1.026 ms 615.6 kcyc   2.2% x1 1026 us/call
  fx_timeline_step        91.38 us  54.8 kcyc   0.2% x1 91.4 us/call
  canvas_clear           100.29 us  60.2 kcyc   0.2% x1 100 us/call
  canvas_buffer_wait       1.00 us    602 cyc   0.0% x1 1.00 us/call
```

Wall min/mean/max was 44.386/47.107/81.767 ms. Render averaged 47.107 ms per frame, including 42.799 ms in `fx_shader_draw`.

This regime leaves 15.393 ms of mean render headroom per display window.

### Window containing the live peak (frames 161–192)

```
frame                    49.092 ms 29.46 Mcyc 100.0%
  pov_preserve_half      137.85 us  82.7 kcyc   0.3% x1 138 us/call
  fx_shader_draw         44.811 ms 26.89 Mcyc  91.3% x1 44811 us/call
    fx_prepare_column    467.63 us 280.6 kcyc   1.0% x144 3.25 us/call
  fx_palette_base_bake   329.38 us 197.6 kcyc   0.7% x1 329 us/call
  fx_prepare_frame        2.450 ms  1.47 Mcyc   5.0% x1 2450 us/call
  fx_advance              1.007 ms 604.1 kcyc   2.1% x1 1007 us/call
    palette_cycle_step   963.29 us 578.0 kcyc   2.0% x1 963 us/call
  fx_timeline_step        96.09 us  57.7 kcyc   0.2% x1 96.1 us/call
  canvas_clear            96.10 us  57.7 kcyc   0.2% x1 96.1 us/call
  canvas_buffer_wait       0.61 us    366 cyc   0.0% x1 0.61 us/call
```

Wall min/mean/max was 46.888/49.091/51.489 ms. Render averaged 49.091 ms per frame, including 44.811 ms in `fx_shader_draw`.

This regime leaves 13.409 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    48.530 ms 29.12 Mcyc 100.0%
  pov_preserve_half      137.60 us  82.6 kcyc   0.3% x1 138 us/call
  fx_shader_draw         44.217 ms 26.53 Mcyc  91.1% x1 44217 us/call
    fx_prepare_column    475.81 us 285.5 kcyc   1.0% x144 3.30 us/call
  fx_palette_base_bake   332.33 us 199.4 kcyc   0.7% x1 332 us/call
  fx_prepare_frame        2.491 ms  1.49 Mcyc   5.1% x1 2491 us/call
  fx_advance              1.004 ms 602.1 kcyc   2.1% x1 1004 us/call
    palette_cycle_step   958.90 us 575.3 kcyc   2.0% x1 959 us/call
  fx_timeline_step        91.87 us  55.1 kcyc   0.2% x1 91.9 us/call
  canvas_clear            96.93 us  58.2 kcyc   0.2% x1 96.9 us/call
  canvas_buffer_wait       0.77 us    462 cyc   0.0% x1 0.77 us/call
```

Wall min/mean/max was 46.043/48.530/50.593 ms. Render averaged 48.529 ms per frame, including 44.217 ms in `fx_shader_draw`.

This regime leaves 13.971 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 161–192, measured over the separate 2,002,189 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1153.3/frame  0.503/1.791/18.251 us  3.30%
  isr_pack           144.2/frame  6.935/7.696/16.470 us  1.77%
  isr_dma_submit     144.2/frame  0.581/0.915/6.631 us  0.21%
```

Inclusive wake time leaves 60.438 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 44.811 ms/frame (91.3% of draw time).
2. `fx_prepare_frame` took 2.450 ms/frame (5.0% of draw time).
3. `fx_advance` took 1.007 ms/frame (2.1% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeMelt` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
