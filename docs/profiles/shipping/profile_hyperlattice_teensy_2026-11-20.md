# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 300-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/hyperlattice_ship.log`. Supersedes `profile_hyperlattice_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | HyperLattice at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 300 s slot, 306 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 306 16` from the captured revision |

Captured 2026-11-20 15:27, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `b55dae00fdf687b1f0e85b62c7fcf77b084ee8a00cd75e9a1d881bd0b85e998c`.

Profile image size `FLASH: code:129528, data:148764, headers:8428 free for files:1744896 / RAM1: variables:313024, code:13832, padding:18936 free for local variables:178496 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `4740005a511fa9822cd203a5d7df0177fd005775ba9c22eff01a6565d3ffa0f0`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 3089–3104. The root counter measured 723481.258 µs, compared with 723,481 µs of draw time (0.36 ppm).

## Frame cadence

**Peak confirmed live render 53.279 ms** at frame 3209 (octet-4d-flight). **0/4,786 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 45.645 ms. Live render averaged 29.859 ms.

Among complete live readouts, the largest window mean was 45.216 ms/frame at frames 3089–3104. The dominant scope averaged 27.775 ms/frame across those readouts.

The shipping ELF assigns 2,400 revolutions (300 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 300 s. Frames 2–4787 are confirmed live. The raw file contains 4,873 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 13 frame(s), 4788–4800, with uncertain membership. Their largest render was 33.948 ms and 0 exceeded the budget. The live denominator lies between 4,786 and 4,799. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 53.279 (9), spilled 🟢 0/4786 (0.00%).

README setup frame 1: 45.645 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    27.047 ms 16.23 Mcyc 100.0%
  pov_preserve_half      133.01 us  79.8 kcyc   0.5% x0.94 142 us/call
  hl_shader_draw         22.331 ms 13.40 Mcyc  82.6% x0.94 23820 us/call
  hl_shader_draw          2.813 ms  1.69 Mcyc  10.4% x0.06 45008 us/call
  palette_cycle_step      1.305 ms 783.3 kcyc   4.8% x1 1305 us/call
  hl_timeline_step        15.32 us   9.2 kcyc   0.1% x1 15.3 us/call
  canvas_clear            96.41 us  57.8 kcyc   0.4% x1 96.4 us/call
  canvas_buffer_wait       1.41 us    848 cyc   0.0% x1 1.41 us/call
```

Wall min/mean/max was 24.966/27.047/45.648 ms. Render averaged 27.046 ms per frame, including 2.813 ms in `hl_shader_draw`.

This regime leaves 35.454 ms of mean render headroom per display window.

### Window containing the live peak (frames 3201–3216)

```
frame                    44.541 ms 26.72 Mcyc 100.0%
  pov_preserve_half      137.30 us  82.4 kcyc   0.3% x1 137 us/call
  hl_shader_draw         42.451 ms 25.47 Mcyc  95.3% x1 42451 us/call
  palette_cycle_step      1.474 ms 884.6 kcyc   3.3% x1 1474 us/call
  hl_timeline_step        15.39 us   9.2 kcyc   0.0% x1 15.4 us/call
  canvas_clear            95.90 us  57.5 kcyc   0.2% x1 95.9 us/call
  canvas_buffer_wait       1.77 us   1.1 kcyc   0.0% x1 1.77 us/call
```

Wall min/mean/max was 38.417/44.540/53.280 ms. Render averaged 44.539 ms per frame, including 42.451 ms in `hl_shader_draw`.

This regime leaves 17.961 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 3089–3104)

```
frame                    45.218 ms 27.13 Mcyc 100.0%
  pov_preserve_half      137.38 us  82.4 kcyc   0.3% x1 137 us/call
  hl_shader_draw         43.161 ms 25.90 Mcyc  95.5% x1 43161 us/call
  palette_cycle_step      1.463 ms 877.8 kcyc   3.2% x1 1463 us/call
  hl_timeline_step         2.55 us   1.5 kcyc   0.0% x1 2.55 us/call
  canvas_clear            94.57 us  56.7 kcyc   0.2% x1 94.6 us/call
  canvas_buffer_wait       1.45 us    868 cyc   0.0% x1 1.45 us/call
```

Wall min/mean/max was 43.625/45.217/46.189 ms. Render averaged 45.216 ms per frame, including 43.161 ms in `hl_shader_draw`.

This regime leaves 17.284 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 3105–3120)

```
frame                    44.732 ms 26.84 Mcyc 100.0%
  pov_preserve_half      138.11 us  82.9 kcyc   0.3% x1 138 us/call
  hl_shader_draw         42.662 ms 25.60 Mcyc  95.4% x1 42662 us/call
  palette_cycle_step      1.462 ms 877.1 kcyc   3.3% x1 1462 us/call
  hl_timeline_step         6.66 us   4.0 kcyc   0.0% x1 6.66 us/call
  canvas_clear            94.90 us  56.9 kcyc   0.2% x1 94.9 us/call
  canvas_buffer_wait       1.43 us    858 cyc   0.0% x1 1.43 us/call
```

Wall min/mean/max was 42.869/44.731/46.558 ms. Render averaged 44.731 ms per frame, including 42.662 ms in `hl_shader_draw`.

This regime leaves 17.769 ms of mean render headroom per display window.

### Final complete in-show window (frames 4769–4784)

```
frame                    33.323 ms 19.99 Mcyc 100.0%
  pov_preserve_half      141.15 us  84.7 kcyc   0.4% x1 141 us/call
  hl_shader_draw         31.282 ms 18.77 Mcyc  93.9% x1 31282 us/call
  palette_cycle_step      1.453 ms 871.7 kcyc   4.4% x1 1453 us/call
  hl_timeline_step         2.57 us   1.5 kcyc   0.0% x1 2.57 us/call
  canvas_clear            87.11 us  52.3 kcyc   0.3% x1 87.1 us/call
  canvas_buffer_wait       1.38 us    830 cyc   0.0% x1 1.38 us/call
```

Wall min/mean/max was 32.730/33.323/33.669 ms. Render averaged 33.322 ms per frame, including 31.282 ms in `hl_shader_draw`.

This regime leaves 29.178 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

Authored departures last 240 frames. Cubic, octet and shell flight interpolate into their wide or close variant; the other departures fade through black. The serial marker lands at interpolation start or fade midpoint. Only complete windows inside the 320-frame dwell qualify.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| octet-4d-flight | held preset | 19 | 3089–3104 | — | — | 43.161 | 45.216 | 53.279 | 16 | 0/559 |
| shell-4d-flight | held preset | 19 | 4593–4608 | — | — | 33.260 | 35.304 | 36.614 | 16 | 0/435 |
| octet-wide-flight | held preset | 19 | 2273–2288 | — | — | 29.362 | 31.340 | 31.947 | 16 | 0/679 |
| hypercube-flight | held preset | 19 | 1185–1200 | — | — | 29.044 | 30.999 | 38.634 | 16 | 0/559 |
| shell-close-flight | held preset | 19 | 4145–4160 | — | — | 29.044 | 31.035 | 34.004 | 16 | 0/679 |
| octet-flight | held preset | 19 | 1809–1824 | — | — | 27.931 | 29.900 | 32.689 | 16 | 0/439 |
| cubic-wide-flight | held preset | 19 | 769–784 | — | — | 25.574 | 27.544 | 29.729 | 16 | 0/679 |
| shell-flight | held preset | 19 | 3553–3568 | — | — | 25.082 | 27.051 | 29.889 | 16 | 0/439 |
| cubic-flight | held preset | 18 | 97–112 | — | — | 24.590 | 26.565 | 27.453 | 16 | 0/318 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 3201–3216, measured over the separate 994,557 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1145.8/frame  0.565/1.593/18.671 us  2.93%
  isr_pack           143.2/frame  4.533/6.130/9.330 us  1.41%
  isr_dma_submit     143.2/frame  0.615/0.916/0.980 us  0.21%
```

Inclusive wake time leaves 60.669 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `hl_shader_draw` took 42.451 ms/frame (95.3% of draw time).
2. `palette_cycle_step` took 1.474 ms/frame (3.3% of draw time).
3. `pov_preserve_half` took 0.137 ms/frame (0.3% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=HyperLattice` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
