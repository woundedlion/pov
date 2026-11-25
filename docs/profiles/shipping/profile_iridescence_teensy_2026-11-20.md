# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/iridescence_ship.log`. Supersedes `profile_iridescence_teensy_2026-11-18.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Iridescence at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Iridescence profile 66 32` from the captured revision |

Captured 2026-11-20 16:17, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `fc5f1d45fb8de9108d8ea8cb5a20510a93353bb5e981b499d66a4e6e45f9f218`.

Profile image size `FLASH: code:112608, data:157428, headers:8492 free for files:1753088 / RAM1: variables:314144, code:26616, padding:6152 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `96fe42fe58e9877f453f4b302b9ec365c1dc5728bf5f4ad159b314a0d6cdeacd`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 289–320. The root counter measured 1452759.822 µs, compared with 1,452,758 µs of draw time (1.25 ppm).

## Frame cadence

**Peak confirmed live render 45.752 ms** at frame 899 (Iridescence). **0/949 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 84.370 ms. Live render averaged 45.286 ms.

Among complete live readouts, the largest window mean was 45.397 ms/frame at frames 289–320. The dominant scope averaged 42.263 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–950 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 10 frame(s), 951–960, with uncertain membership. Their largest render was 45.376 ms and 0 exceeded the budget. The live denominator lies between 949 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 45.752, spilled 🟢 0/949 (0.00%).

README setup frame 1: 84.370 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    46.556 ms 27.93 Mcyc 100.0%
  pov_preserve_half      132.58 us  79.5 kcyc   0.3% x0.97 137 us/call
  fx_shader_draw         43.459 ms 26.08 Mcyc  93.3% x1 43459 us/call
    fx_prepare_column     1.147 ms 688.3 kcyc   2.5% x144 7.97 us/call
  fx_palette_base_bake    75.23 us  45.1 kcyc   0.2% x1 75.2 us/call
  fx_prepare_frame        2.546 ms  1.53 Mcyc   5.5% x1 2546 us/call
  fx_advance              22.86 us  13.7 kcyc   0.0% x1 22.9 us/call
  fx_timeline_step        15.14 us   9.1 kcyc   0.0% x1 15.1 us/call
  canvas_clear            93.92 us  56.4 kcyc   0.2% x1 93.9 us/call
  canvas_buffer_wait       1.78 us   1.1 kcyc   0.0% x1 1.78 us/call
```

Wall min/mean/max was 45.101/46.555/84.370 ms. Render averaged 46.554 ms per frame, including 43.459 ms in `fx_shader_draw`.

This regime leaves 15.946 ms of mean render headroom per display window.

### Window containing the live peak (frames 897–928)

```
frame                    45.355 ms 27.21 Mcyc 100.0%
  pov_preserve_half      139.25 us  83.5 kcyc   0.3% x1 139 us/call
  fx_shader_draw         42.331 ms 25.40 Mcyc  93.3% x1 42331 us/call
    fx_prepare_column     1.164 ms 698.7 kcyc   2.6% x144 8.09 us/call
  fx_palette_base_bake    74.26 us  44.6 kcyc   0.2% x1 74.3 us/call
  fx_prepare_frame        2.481 ms  1.49 Mcyc   5.5% x1 2481 us/call
  fx_advance              21.89 us  13.1 kcyc   0.0% x1 21.9 us/call
  fx_timeline_step        14.90 us   8.9 kcyc   0.0% x1 14.9 us/call
  canvas_clear            89.59 us  53.8 kcyc   0.2% x1 89.6 us/call
  canvas_buffer_wait       1.90 us   1.1 kcyc   0.0% x1 1.90 us/call
```

Wall min/mean/max was 45.139/45.355/45.753 ms. Render averaged 45.353 ms per frame, including 42.331 ms in `fx_shader_draw`.

This regime leaves 17.147 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 897–928, measured over the separate 1,999,915 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1152.0/frame  0.501/1.740/18.951 us  3.20%
  isr_pack           144.0/frame  6.803/7.508/10.165 us  1.73%
  isr_dma_submit     144.0/frame  0.618/0.915/5.726 us  0.21%
```

Inclusive wake time leaves 60.500 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 42.331 ms/frame (93.3% of draw time).
2. `fx_prepare_frame` took 2.481 ms/frame (5.5% of draw time).
3. `fx_prepare_column` took 1.164 ms/frame (2.6% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Iridescence` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
