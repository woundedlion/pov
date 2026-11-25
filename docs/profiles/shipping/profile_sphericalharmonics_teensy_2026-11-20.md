# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 96-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/sphericalharmonics_ship.log`. Supersedes `profile_sphericalharmonics_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | SphericalHarmonics at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 96 s slot, 102 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh SphericalHarmonics profile 102 16` from the captured revision |

Captured 2026-11-20 15:51, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `310c72299ae71f3082cc3d76e7a893b532602aa859b18ec39bed548730b63869`.

Profile image size `FLASH: code:42872, data:145176, headers:8560 free for files:1835008 / RAM1: variables:312896, code:12216, padding:20552 free for local variables:178624 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `3ae4f49c688d815c8388366af844f51ed96f0f228fb728d1406cdec9319f828e`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 1105–1120. The root counter measured 197021.917 µs, compared with 197,021 µs of draw time (4.65 ppm).

## Frame cadence

**Peak confirmed live render 12.368 ms** at frame 1133 (morph toward l=4, m=4). **0/1,519 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 20.591 ms. Live render averaged 11.586 ms.

Among complete live readouts, the largest window mean was 12.312 ms/frame at frames 1105–1120. The dominant scope averaged 11.332 ms/frame across those readouts.

The shipping ELF assigns 768 revolutions (96 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 96 s. Frames 2–1520 are confirmed live. The raw file contains 1,609 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 1521–1536, with uncertain membership. Their largest render was 11.229 ms and 0 exceeded the budget. The live denominator lies between 1,519 and 1,535. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 12.368 (24), spilled 🟢 0/1519 (0.00%).

README setup frame 1: 20.591 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    11.582 ms  6.95 Mcyc 100.0%
  pov_preserve_half      141.36 us  84.8 kcyc   1.2% x0.94 151 us/call
  sh_rasterize           11.309 ms  6.79 Mcyc  97.6% x1 11309 us/call
  sh_timeline_step        29.15 us  17.5 kcyc   0.3% x1 29.1 us/call
  canvas_clear            92.76 us  55.7 kcyc   0.8% x1 92.8 us/call
  canvas_buffer_wait       1.04 us    625 cyc   0.0% x1 1.04 us/call
```

Wall min/mean/max was 10.852/11.581/20.591 ms. Render averaged 11.581 ms per frame, including 11.309 ms in `sh_rasterize`.

This regime leaves 50.919 ms of mean render headroom per display window.

### Window containing the live peak (frames 1121–1136)

```
frame                    12.302 ms  7.38 Mcyc 100.0%
  pov_preserve_half      146.14 us  87.7 kcyc   1.2% x1 146 us/call
  sh_rasterize           12.040 ms  7.22 Mcyc  97.9% x1 12040 us/call
  sh_timeline_step        23.30 us  14.0 kcyc   0.2% x1 23.3 us/call
  canvas_clear            83.66 us  50.2 kcyc   0.7% x1 83.7 us/call
  canvas_buffer_wait       1.25 us    748 cyc   0.0% x1 1.25 us/call
```

Wall min/mean/max was 12.208/12.302/12.369 ms. Render averaged 12.301 ms per frame, including 12.040 ms in `sh_rasterize`.

This regime leaves 50.199 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 673–688)

```
frame                    12.300 ms  7.38 Mcyc 100.0%
  pov_preserve_half      144.12 us  86.5 kcyc   1.2% x1 144 us/call
  sh_rasterize           12.044 ms  7.23 Mcyc  97.9% x1 12044 us/call
  sh_timeline_step        19.91 us  11.9 kcyc   0.2% x1 19.9 us/call
  canvas_clear            83.00 us  49.8 kcyc   0.7% x1 83.0 us/call
  canvas_buffer_wait       1.48 us    888 cyc   0.0% x1 1.48 us/call
```

Wall min/mean/max was 12.255/12.300/12.344 ms. Render averaged 12.299 ms per frame, including 12.044 ms in `sh_rasterize`.

This regime leaves 50.201 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 1105–1120)

```
frame                    12.314 ms  7.39 Mcyc 100.0%
  pov_preserve_half      143.48 us  86.1 kcyc   1.2% x1 143 us/call
  sh_rasterize           12.059 ms  7.24 Mcyc  97.9% x1 12059 us/call
  sh_timeline_step        19.81 us  11.9 kcyc   0.2% x1 19.8 us/call
  canvas_clear            82.91 us  49.7 kcyc   0.7% x1 82.9 us/call
  canvas_buffer_wait       1.48 us    888 cyc   0.0% x1 1.48 us/call
```

Wall min/mean/max was 12.283/12.313/12.359 ms. Render averaged 12.312 ms per frame, including 12.059 ms in `sh_rasterize`.

This regime leaves 50.188 ms of mean render headroom per display window.

### Final complete in-show window (frames 1505–1520)

```
frame                    11.019 ms  6.61 Mcyc 100.0%
  pov_preserve_half      148.82 us  89.3 kcyc   1.4% x1 149 us/call
  sh_rasterize           10.754 ms  6.45 Mcyc  97.6% x1 10754 us/call
  sh_timeline_step        23.50 us  14.1 kcyc   0.2% x1 23.5 us/call
  canvas_clear            83.42 us  50.1 kcyc   0.8% x1 83.4 us/call
  canvas_buffer_wait       1.25 us    748 cyc   0.0% x1 1.25 us/call
```

Wall min/mean/max was 10.876/11.018/11.181 ms. Render averaged 11.017 ms per frame, including 10.754 ms in `sh_rasterize`.

This regime leaves 51.483 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

The mode marker names the harmonic target. Each representative window lies in the later half of a continuous 64-frame morph.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| l=4, m=-3 | late morph toward target | 1 | 673–688 | — | — | 12.044 | 12.299 | 12.343 | 16 | 0/64 |
| l=4, m=4 | late morph toward target | 1 | 1121–1136 | — | — | 12.040 | 12.301 | 12.368 | 16 | 0/64 |
| l=4, m=-2 | late morph toward target | 1 | 737–752 | — | — | 11.882 | 12.144 | 12.229 | 16 | 0/64 |
| l=4, m=-4 | late morph toward target | 1 | 609–624 | — | — | 11.871 | 12.130 | 12.251 | 16 | 0/64 |
| l=4, m=3 | late morph toward target | 1 | 1057–1072 | — | — | 11.859 | 12.115 | 12.176 | 16 | 0/64 |
| l=4, m=-1 | late morph toward target | 1 | 801–816 | — | — | 11.688 | 11.943 | 12.000 | 16 | 0/64 |
| l=4, m=2 | late morph toward target | 1 | 993–1008 | — | — | 11.685 | 11.939 | 12.042 | 16 | 0/64 |
| l=3, m=-2 | late morph toward target | 1 | 225–240 | — | — | 11.571 | 11.831 | 11.898 | 16 | 0/64 |
| l=3, m=3 | late morph toward target | 1 | 545–560 | — | — | 11.544 | 11.806 | 11.878 | 16 | 0/64 |
| l=3, m=-3 | late morph toward target | 1 | 161–176 | — | — | 11.486 | 11.752 | 11.799 | 16 | 0/64 |
| l=1, m=-1 | late morph toward target | 1 | 1185–1200 | — | — | 11.442 | 11.701 | 12.197 | 16 | 0/64 |
| l=4, m=0 | late morph toward target | 1 | 865–880 | — | — | 11.416 | 11.671 | 11.996 | 16 | 0/64 |
| l=3, m=2 | late morph toward target | 1 | 481–496 | — | — | 11.411 | 11.668 | 11.784 | 16 | 0/64 |
| l=3, m=-1 | late morph toward target | 1 | 289–304 | — | — | 11.364 | 11.620 | 11.751 | 16 | 0/64 |
| l=4, m=1 | late morph toward target | 1 | 929–944 | — | — | 11.354 | 11.609 | 11.721 | 16 | 0/64 |
| l=3, m=0 | late morph toward target | 1 | 353–368 | — | — | 11.224 | 11.477 | 11.586 | 16 | 0/64 |
| l=3, m=1 | late morph toward target | 1 | 417–432 | — | — | 11.142 | 11.400 | 11.521 | 16 | 0/64 |
| l=2, m=2 | late morph toward target | 1 | 97–112 | — | — | 11.095 | 11.362 | 11.456 | 16 | 0/64 |
| l=2, m=-1 | late morph toward target | 1 | 1441–1456 | — | — | 11.078 | 11.343 | 11.416 | 16 | 0/64 |
| l=2, m=-2 | late morph toward target | 1 | 1377–1392 | — | — | 11.023 | 11.289 | 11.371 | 16 | 0/64 |
| l=2, m=1 | late morph toward target | 1 | 33–48 | — | — | 10.827 | 11.084 | 11.234 | 16 | 0/62 |
| l=1, m=1 | late morph toward target | 1 | 1313–1328 | — | — | 10.575 | 10.838 | 10.919 | 16 | 0/64 |
| l=1, m=0 | late morph toward target | 1 | 1249–1264 | — | — | 10.571 | 10.833 | 11.543 | 16 | 0/64 |
| l=2, m=0 | late morph toward target | 0 | no qualified window | — | — | — | — | 11.225 | 16 | 0/49 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 1121–1136, measured over the separate 999,865 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.9/frame  0.446/1.670/12.831 us  3.07%
  isr_pack           144.0/frame  6.853/7.478/10.150 us  1.72%
  isr_dma_submit     144.0/frame  0.645/0.918/0.978 us  0.21%
```

Inclusive wake time leaves 60.581 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `sh_rasterize` took 12.040 ms/frame (97.9% of draw time).
2. `pov_preserve_half` took 0.146 ms/frame (1.2% of draw time).
3. `canvas_clear` took 0.084 ms/frame (0.7% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/shader.h` and the pullback runtime/palette kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=SphericalHarmonics` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
