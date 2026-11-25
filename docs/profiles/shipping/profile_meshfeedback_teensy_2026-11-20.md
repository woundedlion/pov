# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 181-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/meshfeedback_ship.log`. Supersedes `profile_meshfeedback_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | MeshFeedback at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 181 s slot, 187 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh MeshFeedback profile 187 16` from the captured revision |

Captured 2026-11-20 15:40, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `1e346f10297d3b6471eb4eb04dd019158da5b0d59a2498bf6313cc85a760a62c`.

Profile image size `FLASH: code:194856, data:190408, headers:8976 free for files:1637376 / RAM1: variables:314176, code:45464, padding:20072 free for local variables:144576 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `8a9f5668865e8092997d3417f5f6c86fd1e5584455b142e53dfc74bc1eba886f`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 513–528. The root counter measured 830616.668 µs, compared with 830,616 µs of draw time (0.80 ppm).

## Frame cadence

**Peak confirmed live render 55.879 ms** at frame 2651 (WigglingWormhole). **0/2,886 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 6.551 ms. Live render averaged 47.211 ms.

Among complete live readouts, the largest window mean was 53.279 ms/frame at frames 2049–2064. The dominant scope averaged 38.586 ms/frame across those readouts.

The shipping ELF assigns 1,448 revolutions (181 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 181 s. Frames 2–2887 are confirmed live. The raw file contains 2,968 frames, including 76 confirmed after the boundary.

Boundary resolution leaves 5 frame(s), 2888–2892, with uncertain membership. Their largest render was 52.733 ms and 0 exceeded the budget. The live denominator lies between 2,886 and 2,891. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers the full canvas (41,472 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 55.879 (12), spilled 🟢 0/2886 (0.00%).

README setup frame 1: 6.551 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    41.424 ms 24.85 Mcyc 100.0%
  feedback_forward        1.566 ms 939.9 kcyc   3.8% x0.94 1671 us/call
  mf_mesh_draw            6.783 ms  4.07 Mcyc  16.4% x1 6783 us/call
  mf_feedback_flush      32.904 ms 19.74 Mcyc  79.4% x1 32904 us/call
    feedback_composite   24.176 ms 14.51 Mcyc  58.4% x0.94 25788 us/call
    feedback_populate     7.245 ms  4.35 Mcyc  17.5% x0.94 7728 us/call
    feedback_litscan      79.04 us  47.4 kcyc   0.2% x1 79.0 us/call
  mf_timeline_step        99.62 us  59.8 kcyc   0.2% x1 99.6 us/call
    animation_noise_clock   0.18 us    107 cyc   0.0% x1 0.18 us/call
  mf_apply_params          2.91 us   1.7 kcyc   0.0% x1 2.91 us/call
  canvas_buffer_wait       2.00 us   1.2 kcyc   0.0% x1 2.00 us/call
```

Wall min/mean/max was 6.554/41.424/46.930 ms. Render averaged 41.422 ms per frame, including 32.904 ms in `mf_feedback_flush`.

This regime leaves 21.078 ms of mean render headroom per display window.

### Window containing the live peak (frames 2641–2656)

```
frame                    44.750 ms 26.85 Mcyc 100.0%
  feedback_forward       783.42 us 470.0 kcyc   1.8% x1 783 us/call
  mf_mesh_draw            9.479 ms  5.69 Mcyc  21.2% x1 9479 us/call
  mf_feedback_flush      34.275 ms 20.56 Mcyc  76.6% x1 34275 us/call
    feedback_composite   33.253 ms 19.95 Mcyc  74.3% x1 33253 us/call
    feedback_populate    493.36 us 296.0 kcyc   1.1% x1 493 us/call
    feedback_litscan      12.03 us   7.2 kcyc   0.0% x1 12.0 us/call
  mf_timeline_step        91.54 us  54.9 kcyc   0.2% x1 91.5 us/call
  mf_apply_params         48.64 us  29.2 kcyc   0.1% x1 48.6 us/call
  canvas_buffer_wait       8.06 us   4.8 kcyc   0.0% x1 8.06 us/call
```

Wall min/mean/max was 42.474/44.749/55.881 ms. Render averaged 44.742 ms per frame, including 34.275 ms in `mf_feedback_flush`.

This regime leaves 17.758 ms of mean render headroom per display window.

### Most expensive qualified preset window (frames 513–528)

```
frame                    51.914 ms 31.15 Mcyc 100.0%
  feedback_forward        1.100 ms 659.8 kcyc   2.1% x1 1100 us/call
  mf_mesh_draw            7.756 ms  4.65 Mcyc  14.9% x1 7756 us/call
  mf_feedback_flush      42.901 ms 25.74 Mcyc  82.6% x1 42901 us/call
    feedback_composite   33.574 ms 20.14 Mcyc  64.7% x1 33574 us/call
    feedback_populate     7.804 ms  4.68 Mcyc  15.0% x1 7804 us/call
    feedback_litscan      14.53 us   8.7 kcyc   0.0% x1 14.5 us/call
  mf_timeline_step        90.14 us  54.1 kcyc   0.2% x1 90.1 us/call
    animation_noise_clock   0.18 us    106 cyc   0.0% x1 0.18 us/call
  mf_apply_params          2.61 us   1.6 kcyc   0.0% x1 2.61 us/call
  canvas_buffer_wait       1.82 us   1.1 kcyc   0.0% x1 1.82 us/call
```

Wall min/mean/max was 51.361/51.913/52.458 ms. Render averaged 51.912 ms per frame, including 42.901 ms in `mf_feedback_flush`.

This regime leaves 10.588 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 1201–1216)

```
frame                    51.460 ms 30.88 Mcyc 100.0%
  feedback_forward        1.336 ms 801.8 kcyc   2.6% x1 1336 us/call
  mf_mesh_draw            7.954 ms  4.77 Mcyc  15.5% x1 7954 us/call
  mf_feedback_flush      41.979 ms 25.19 Mcyc  81.6% x1 41979 us/call
    feedback_composite   32.709 ms 19.63 Mcyc  63.6% x1 32709 us/call
    feedback_populate     7.758 ms  4.65 Mcyc  15.1% x1 7758 us/call
    feedback_litscan      12.57 us   7.5 kcyc   0.0% x1 12.6 us/call
  mf_timeline_step        93.19 us  55.9 kcyc   0.2% x1 93.2 us/call
    animation_noise_clock   0.18 us    106 cyc   0.0% x1 0.18 us/call
  mf_apply_params         22.66 us  13.6 kcyc   0.0% x1 22.7 us/call
  canvas_buffer_wait       7.47 us   4.5 kcyc   0.0% x1 7.47 us/call
```

Wall min/mean/max was 49.736/51.460/53.680 ms. Render averaged 51.453 ms per frame, including 41.979 ms in `mf_feedback_flush`.

This regime leaves 11.047 ms of mean render headroom per display window.

### Final complete in-show window (frames 2865–2880)

```
frame                    47.848 ms 28.71 Mcyc 100.0%
  feedback_forward        1.036 ms 621.3 kcyc   2.2% x1 1036 us/call
  mf_mesh_draw           11.867 ms  7.12 Mcyc  24.8% x1 11867 us/call
  mf_feedback_flush      34.783 ms 20.87 Mcyc  72.7% x1 34783 us/call
    feedback_composite   34.266 ms 20.56 Mcyc  71.6% x1 34266 us/call
    feedback_populate     29.49 us  17.7 kcyc   0.1% x1 29.5 us/call
    feedback_litscan      22.74 us  13.6 kcyc   0.0% x1 22.7 us/call
  mf_timeline_step        95.59 us  57.4 kcyc   0.2% x1 95.6 us/call
  mf_apply_params          4.31 us   2.6 kcyc   0.0% x1 4.31 us/call
  canvas_buffer_wait       2.04 us   1.2 kcyc   0.0% x1 2.04 us/call
```

Wall min/mean/max was 46.845/47.848/49.304 ms. Render averaged 47.846 ms per frame, including 34.783 ms in `mf_feedback_flush`.

This regime leaves 14.654 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| EnergeticFire | held preset | 14 | 513–528 | — | — | 42.901 | 51.912 | 53.098 | 16 | 0/241 |
| Miasma | held preset | 14 | 2129–2144 | — | — | 42.772 | 52.961 | 53.937 | 16 | 0/241 |
| SlowDust | held preset | 14 | 1137–1152 | — | — | 42.612 | 51.441 | 52.194 | 16 | 0/241 |
| Smoke | held preset | 14 | 737–752 | — | — | 41.787 | 45.349 | 46.179 | 16 | 0/241 |
| SlowFire | held preset | 14 | 257–272 | — | — | 41.656 | 47.511 | 48.470 | 16 | 0/241 |
| WavyTrails | held preset | 14 | 1297–1312 | — | — | 41.333 | 51.069 | 53.664 | 16 | 0/241 |
| MeltingLo | held preset | 14 | 1809–1824 | — | — | 40.978 | 51.140 | 53.445 | 16 | 0/241 |
| MeltingHi | held preset | 14 | 1569–1584 | — | — | 38.731 | 44.324 | 48.720 | 16 | 0/241 |
| ArcingLightning | held preset | 13 | 17–32 | — | — | 37.422 | 46.767 | 48.441 | 16 | 0/239 |
| LooseWormhole | held preset | 14 | 2369–2384 | — | — | 35.045 | 42.508 | 50.958 | 16 | 0/241 |
| WigglingWormhole | held preset | 14 | 2865–2880 | — | — | 34.783 | 47.846 | 55.879 | 16 | 0/237 |
| TightWormhole | held preset | 14 | 2513–2528 | — | — | 34.651 | 44.853 | 52.754 | 16 | 0/241 |

### Per-pixel figures

The representative window has no complete blend counter for its 41,472 output pixels.

## Column-ISR / DMA marshaling cost

Frames 2641–2656, measured over the separate 1,003,868 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1156.5/frame  0.478/1.801/16.338 us  3.32%
  isr_pack           144.6/frame  6.720/7.941/10.251 us  1.83%
  isr_dma_submit     144.6/frame  0.755/0.918/5.765 us  0.21%
```

Inclusive wake time leaves 60.425 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `mf_feedback_flush` took 34.275 ms/frame (76.6% of draw time).
2. `feedback_composite` took 33.253 ms/frame (74.3% of draw time).
3. `mf_mesh_draw` took 9.479 ms/frame (21.2% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=MeshFeedback` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
