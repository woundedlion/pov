# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/hopffibration_ship.log`. Supersedes `profile_hopffibration_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | HopfFibration at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HopfFibration profile 66 32` from the captured revision |

Captured 2026-11-20 16:15, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `49122614846a75246aa84a0e558898df5dec82c08f4d3e43b9f9257a9b0fc956`.

Profile image size `FLASH: code:130600, data:145680, headers:8392 free for files:1746944 / RAM1: variables:312992, code:29192, padding:3576 free for local variables:178528 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `e87e7679b4f350e5d5ddc08fa76f087d3925f03b02dc0ff774196f9ae361dbab`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 193–224. The root counter measured 1497883.832 µs, compared with 1,497,883 µs of draw time (0.56 ppm).

## Frame cadence

**Peak confirmed live render 55.490 ms** at frame 822 (HopfFibration). **0/945 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 0.672 ms. Live render averaged 26.679 ms.

Among complete live readouts, the largest window mean was 46.807 ms/frame at frames 193–224. The dominant scope averaged 26.687 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–946 are confirmed live. The raw file contains 1,033 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 14 frame(s), 947–960, with uncertain membership. Their largest render was 42.154 ms and 0 exceeded the budget. The live denominator lies between 945 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 55.490, spilled 🟢 0/945 (0.00%).

README setup frame 1: 0.672 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    11.198 ms  6.72 Mcyc 100.0%
  pov_preserve_half      141.50 us  84.9 kcyc   1.3% x0.97 146 us/call
  hf_render_trails       10.731 ms  6.44 Mcyc  95.8% x1 10731 us/call
  hf_project_record      151.96 us  91.2 kcyc   1.4% x1 152 us/call
  hf_advance_tumble       12.00 us   7.2 kcyc   0.1% x1 12.0 us/call
  hf_timeline_step        52.63 us  31.6 kcyc   0.5% x1 52.6 us/call
  canvas_clear            92.40 us  55.4 kcyc   0.8% x1 92.4 us/call
  canvas_buffer_wait       0.97 us    585 cyc   0.0% x1 0.97 us/call
```

Wall min/mean/max was 0.672/11.198/24.343 ms. Render averaged 11.197 ms per frame, including 10.731 ms in `hf_render_trails`.

This regime leaves 51.303 ms of mean render headroom per display window.

### Window containing the live peak (frames 801–832)

```
frame                    40.042 ms 24.03 Mcyc 100.0%
  pov_preserve_half      138.46 us  83.1 kcyc   0.3% x1 138 us/call
  hf_render_trails       39.522 ms 23.71 Mcyc  98.7% x1 39522 us/call
  hf_project_record      186.46 us 111.9 kcyc   0.5% x1 186 us/call
  hf_advance_tumble       12.29 us   7.4 kcyc   0.0% x1 12.3 us/call
  hf_timeline_step        73.17 us  43.9 kcyc   0.2% x1 73.2 us/call
  canvas_clear            92.28 us  55.4 kcyc   0.2% x1 92.3 us/call
  canvas_buffer_wait       1.52 us    915 cyc   0.0% x1 1.52 us/call
```

Wall min/mean/max was 12.468/40.041/55.490 ms. Render averaged 40.040 ms per frame, including 39.522 ms in `hf_render_trails`.

This regime leaves 22.460 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    19.343 ms 11.61 Mcyc 100.0%
  pov_preserve_half      143.15 us  85.9 kcyc   0.7% x1 143 us/call
  hf_render_trails       18.846 ms 11.31 Mcyc  97.4% x1 18846 us/call
  hf_project_record      178.68 us 107.2 kcyc   0.9% x1 179 us/call
  hf_advance_tumble       12.26 us   7.4 kcyc   0.1% x1 12.3 us/call
  hf_timeline_step        55.82 us  33.5 kcyc   0.3% x1 55.8 us/call
  canvas_clear            90.37 us  54.2 kcyc   0.5% x1 90.4 us/call
  canvas_buffer_wait       1.52 us    914 cyc   0.0% x1 1.52 us/call
```

Wall min/mean/max was 7.824/19.343/35.157 ms. Render averaged 19.342 ms per frame, including 18.846 ms in `hf_render_trails`.

This regime leaves 43.158 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 801–832, measured over the separate 2,009,236 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1157.4/frame  0.466/1.605/23.240 us  2.96%
  isr_pack           144.7/frame  4.611/6.396/10.025 us  1.47%
  isr_dma_submit     144.7/frame  0.666/0.926/9.615 us  0.21%
```

Inclusive wake time leaves 60.650 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `hf_render_trails` took 39.522 ms/frame (98.7% of draw time).
2. `hf_project_record` took 0.186 ms/frame (0.5% of draw time).
3. `pov_preserve_half` took 0.138 ms/frame (0.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared scan-shapes, spherical-math and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=HopfFibration` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
