# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 80-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/mindsplatter_ship.log`. Supersedes `profile_mindsplatter_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | MindSplatter at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 80 s slot, 86 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh MindSplatter profile 86 16` from the captured revision |

Captured 2026-11-20 16:09, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `4f21691480a59e3bf302b3085bda65de65ebb4d66cdab90e3f37183e43282fc0`.

Profile image size `FLASH: code:160448, data:201932, headers:8284 free for files:1660952 / RAM1: variables:314208, code:41272, padding:24264 free for local variables:144544 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `3ae34316468400ec448c79cf5b3f000c9ccd3d7ac8d4e9c3f86806f1c7498c26`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 993–1008. The root counter measured 656251.533 µs, compared with 656,250 µs of draw time (2.34 ppm).

## Frame cadence

**Peak confirmed live render 45.451 ms** at frame 1003 (6 → 7). **0/1,263 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 0.423 ms. Live render averaged 21.631 ms.

Among complete live readouts, the largest window mean was 41.014 ms/frame at frames 993–1008. The dominant scope averaged 18.921 ms/frame across those readouts.

The shipping ELF assigns 640 revolutions (80 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 80 s. Frames 2–1264 are confirmed live. The raw file contains 1,353 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 1265–1280, with uncertain membership. Their largest render was 27.015 ms and 0 exceeded the budget. The live denominator lies between 1,263 and 1,279. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 45.451 (8), spilled 🟢 0/1263 (0.00%).

README setup frame 1: 0.423 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    933.00 us 559.8 kcyc 100.0%
  pov_preserve_half      138.45 us  83.1 kcyc  14.8% x0.94 148 us/call
  msp_draw_particles     513.47 us 308.1 kcyc  55.0% x1 513 us/call
    msp_particle_scan      1.41 us    845 cyc   0.2% x0.06 22.5 us/call
  msp_particle_step      128.22 us  76.9 kcyc  13.7% x1 128 us/call
  msp_timeline_step       45.13 us  27.1 kcyc   4.8% x1 45.1 us/call
  canvas_clear            96.25 us  57.8 kcyc  10.3% x1 96.3 us/call
  canvas_buffer_wait       0.38 us    229 cyc   0.0% x1 0.38 us/call
```

Wall min/mean/max was 0.424/0.933/1.375 ms. Render averaged 0.933 ms per frame, including 0.513 ms in `msp_draw_particles`.

This regime leaves 61.567 ms of mean render headroom per display window.

### Window containing the live peak (frames 993–1008)

```
frame                    41.016 ms 24.61 Mcyc 100.0%
  pov_preserve_half      138.10 us  82.9 kcyc   0.3% x1 138 us/call
  msp_draw_particles     37.573 ms 22.54 Mcyc  91.6% x1 37573 us/call
  msp_particle_step       3.059 ms  1.84 Mcyc   7.5% x1 3059 us/call
  msp_timeline_step      130.93 us  78.6 kcyc   0.3% x1 131 us/call
  canvas_clear            95.41 us  57.2 kcyc   0.2% x1 95.4 us/call
  canvas_buffer_wait       1.39 us    832 cyc   0.0% x1 1.39 us/call
```

Wall min/mean/max was 37.482/41.015/45.452 ms. Render averaged 41.014 ms per frame, including 37.573 ms in `msp_draw_particles`.

This regime leaves 21.486 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 1009–1024)

```
frame                    33.933 ms 20.36 Mcyc 100.0%
  pov_preserve_half      136.03 us  81.6 kcyc   0.4% x1 136 us/call
  msp_draw_particles     30.781 ms 18.47 Mcyc  90.7% x1 30781 us/call
  msp_particle_step       2.797 ms  1.68 Mcyc   8.2% x1 2797 us/call
  msp_timeline_step      104.06 us  62.4 kcyc   0.3% x1 104 us/call
  canvas_clear            95.46 us  57.3 kcyc   0.3% x1 95.5 us/call
  canvas_buffer_wait       1.46 us    874 cyc   0.0% x1 1.46 us/call
```

Wall min/mean/max was 27.019/33.932/39.759 ms. Render averaged 33.931 ms per frame, including 30.781 ms in `msp_draw_particles`.

This regime leaves 28.569 ms of mean render headroom per display window.

### Final complete in-show window (frames 1249–1264)

```
frame                    23.530 ms 14.12 Mcyc 100.0%
  pov_preserve_half      139.24 us  83.5 kcyc   0.6% x1 139 us/call
  msp_draw_particles     21.021 ms 12.61 Mcyc  89.3% x1 21021 us/call
  msp_particle_step       2.155 ms  1.29 Mcyc   9.2% x1 2155 us/call
  msp_timeline_step      101.54 us  60.9 kcyc   0.4% x1 102 us/call
  canvas_clear            95.24 us  57.1 kcyc   0.4% x1 95.2 us/call
  canvas_buffer_wait       1.39 us    832 cyc   0.0% x1 1.39 us/call
```

Wall min/mean/max was 21.138/23.530/25.735 ms. Render averaged 23.529 ms per frame, including 21.021 ms in `msp_draw_particles`.

This regime leaves 38.971 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| 7 | held preset | 6 | 1009–1024 | — | — | 30.781 | 33.931 | 45.451 | 16 | 0/159 |
| 6 | held preset | 6 | 913–928 | — | — | 29.902 | 37.993 | 42.990 | 16 | 0/159 |
| 4 | held preset | 6 | 529–544 | — | — | 26.222 | 29.626 | 44.171 | 16 | 0/159 |
| 2 | held preset | 6 | 209–224 | — | — | 25.458 | 28.113 | 43.447 | 16 | 0/159 |
| 8 | held preset | 5 | 1233–1248 | — | — | 21.519 | 24.057 | 27.086 | 16 | 0/151 |
| 3 | held preset | 6 | 369–384 | — | — | 19.650 | 22.136 | 25.878 | 16 | 0/159 |
| 1 | held preset | 8 | 129–144 | — | — | 14.392 | 16.751 | 19.041 | 16 | 0/158 |
| 5 | held preset | 6 | 689–704 | — | — | 13.770 | 15.452 | 26.999 | 16 | 0/159 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 993–1008, measured over the separate 1,002,099 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1154.4/frame  0.590/1.663/15.551 us  3.06%
  isr_pack           144.3/frame  4.791/6.638/9.365 us  1.52%
  isr_dma_submit     144.3/frame  0.673/0.908/1.186 us  0.20%
```

Inclusive wake time leaves 60.588 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `msp_draw_particles` took 37.573 ms/frame (91.6% of draw time).
2. `msp_particle_step` took 3.059 ms/frame (7.5% of draw time).
3. `pov_preserve_half` took 0.138 ms/frame (0.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include the shared Plot raster, column-walk and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=MindSplatter` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
