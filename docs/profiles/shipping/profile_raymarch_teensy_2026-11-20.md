# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 30-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/raymarch_ship.log`. Supersedes `profile_raymarch_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | Raymarch at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 30 s slot, 36 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Raymarch profile 36 32` from the captured revision |

Captured 2026-11-20 16:05, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `2e493b88f576e9935eb7445f2a11d5c34bbad722c9853c89fa7be9285492faf1`.

Profile image size `FLASH: code:105760, data:186396, headers:8900 free for files:1730560 / RAM1: variables:313920, code:25976, padding:6792 free for local variables:177600 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `86299bec1876f5ce22ca1b8d2ea2c4997994799b04c9632108ab3953229bc728`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 321–352. The root counter measured 1532900.093 µs, compared with 1,532,901 µs of draw time (0.59 ppm).

## Frame cadence

**Peak confirmed live render 52.092 ms** at frame 10 (Raymarch). **0/469 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 93.629 ms. Live render averaged 46.912 ms.

Among complete live readouts, the largest window mean was 47.902 ms/frame at frames 321–352. The dominant scope averaged 44.075 ms/frame across those readouts.

The shipping ELF assigns 240 revolutions (30 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 30 s. Frames 2–470 are confirmed live. The raw file contains 552 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 10 frame(s), 471–480, with uncertain membership. Their largest render was 49.564 ms and 0 exceeded the budget. The live denominator lies between 469 and 479. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 52.092, spilled 🟢 0/469 (0.00%).

README setup frame 1: 93.629 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    48.734 ms 29.24 Mcyc 100.0%
  pov_preserve_half      128.82 us  77.3 kcyc   0.3% x0.97 133 us/call
  rm_shader_draw         45.890 ms 27.53 Mcyc  94.2% x1 45890 us/call
  rm_timeline_step       350.71 us 210.4 kcyc   0.7% x1 351 us/call
  canvas_clear            89.95 us  54.0 kcyc   0.2% x1 90.0 us/call
  canvas_buffer_wait       1.77 us   1.1 kcyc   0.0% x1 1.77 us/call
```

Wall min/mean/max was 43.261/48.734/93.635 ms. Render averaged 48.733 ms per frame, including 45.890 ms in `rm_shader_draw`.

This regime leaves 13.767 ms of mean render headroom per display window.

### Final complete in-show window (frames 417–448)

```
frame                    47.310 ms 28.39 Mcyc 100.0%
  pov_preserve_half      133.77 us  80.3 kcyc   0.3% x1 134 us/call
  rm_shader_draw         44.474 ms 26.68 Mcyc  94.0% x1 44474 us/call
  rm_timeline_step       346.77 us 208.1 kcyc   0.7% x1 347 us/call
  canvas_clear            85.44 us  51.3 kcyc   0.2% x1 85.4 us/call
  canvas_buffer_wait       1.73 us   1.0 kcyc   0.0% x1 1.73 us/call
```

Wall min/mean/max was 44.451/47.309/49.151 ms. Render averaged 47.308 ms per frame, including 44.474 ms in `rm_shader_draw`.

This regime leaves 15.192 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 1–32, measured over the separate 2,417,117 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1392.3/frame  0.471/1.391/12.718 us  2.56%
  isr_pack           138.2/frame  4.536/6.273/9.683 us  1.14%
  isr_dma_submit     138.2/frame  0.571/0.913/1.193 us  0.16%
```

Inclusive wake time leaves 60.900 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `rm_shader_draw` took 45.890 ms/frame (94.2% of draw time).
2. `rm_timeline_step` took 0.351 ms/frame (0.7% of draw time).
3. `pov_preserve_half` took 0.129 ms/frame (0.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=Raymarch` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
