# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 306-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/islamicstars_ship.log`. Supersedes `profile_islamicstars_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3 (serial 6362410), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | IslamicStars at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 306 s slot, 312 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM3 bash tools/profile_one.sh IslamicStars profile 312 16` from the captured revision |

Captured 2026-11-20 15:00, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `55f1d04f57945d5c692c898f2622cdf66c0795445d7d35aa203cecb0b337754b`.

Profile image size `FLASH: code:184928, data:196032, headers:9184 free for files:1641472 / RAM1: variables:314688, code:30728, padding:2040 free for local variables:176832 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `2150d8e96f8d0e6aa1d3f0deb71c5286d30f2ed8a79c4b29174461a129fcd456`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 2049–2064. The root counter measured 906072.893 µs, compared with 906,070 µs of draw time (3.19 ppm).

## Frame cadence

**Peak confirmed live render 58.986 ms** at frame 2056 (dodecahedron_hk35_ambo_hk62_ambo_relax_hk42). **0/4,880 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 16.611 ms. Live render averaged 24.907 ms.

Among complete live readouts, the largest window mean was 56.629 ms/frame at frames 2049–2064. The dominant scope averaged 24.661 ms/frame across those readouts.

The shipping ELF assigns 2,448 revolutions (306 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 306 s. Frames 2–4881 are confirmed live. The raw file contains 4,969 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 15 frame(s), 4882–4896, with uncertain membership. Their largest render was 27.467 ms and 0 exceeded the budget. The live denominator lies between 4,880 and 4,895. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 58.986 (16), spilled 🟢 0/4880 (0.00%).

README setup frame 1: 16.611 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    10.546 ms  6.33 Mcyc 100.0%
  pov_preserve_half      138.31 us  83.0 kcyc   1.3% x0.94 148 us/call
  is_timeline_step       10.290 ms  6.17 Mcyc  97.6% x1 10290 us/call
    is_draw_shape        10.114 ms  6.07 Mcyc  95.9% x1 10114 us/call
      is_mesh_scan        9.827 ms  5.90 Mcyc  93.2% x1 9827 us/call
        scan_mesh_raster  9.416 ms  5.65 Mcyc  89.3% x7.94 1186 us/call
        scan_face_setup  217.39 us 130.4 kcyc   2.1% x7.94 27.4 us/call
      is_face_offsets     18.30 us  11.0 kcyc   0.2% x1 18.3 us/call
      is_mesh_transform   20.14 us  12.1 kcyc   0.2% x1 20.1 us/call
      mesh_visibility     23.84 us  14.3 kcyc   0.2% x1 23.8 us/call
  is_ripple_prepare        6.44 us   3.9 kcyc   0.1% x1 6.44 us/call
  canvas_clear            94.64 us  56.8 kcyc   0.9% x1 94.6 us/call
  canvas_buffer_wait       0.13 us     77 cyc   0.0% x1 0.13 us/call
```

Wall min/mean/max was 8.727/10.546/16.611 ms. Render averaged 10.546 ms per frame, including 10.290 ms in `is_timeline_step`.

This regime leaves 51.954 ms of mean render headroom per display window.

### Window containing the live peak (frames 2049–2064)

```
scan_face_setup          363.37 us 218.0 kcyc   0.6% x17.88 20.3 us/call
mesh_copy_capture        192.24 us 115.3 kcyc   0.3% x22.19 8.66 us/call
mesh_capture              2.617 ms  1.57 Mcyc   4.6% x5 523 us/call
mesh_classify             1.493 ms 896.1 kcyc   2.6% x1 1493 us/call
frame                    56.630 ms 33.98 Mcyc 100.0%
  pov_preserve_half      141.01 us  84.6 kcyc   0.2% x1 141 us/call
  is_timeline_step       56.324 ms 33.79 Mcyc  99.5% x1 56324 us/call
    is_draw_shape        56.146 ms 33.69 Mcyc  99.1% x1 56146 us/call
      is_mesh_scan       46.808 ms 28.08 Mcyc  82.7% x1 46808 us/call
        scan_face_sweep   5.504 ms  3.30 Mcyc   9.7% x1 5504 us/call
          scan_sweep_column  2.303 ms  1.38 Mcyc   4.1% x143.25 16.1 us/call
        mesh_sweep_paths 109.87 us  65.9 kcyc   0.2% x1 110 us/call
        scan_mesh_raster 27.429 ms 16.46 Mcyc  48.4% x387.88 70.7 us/call
        scan_face_setup   5.563 ms  3.34 Mcyc   9.8% x426.88 13.0 us/call
      is_face_offsets     44.77 us  26.9 kcyc   0.1% x1 44.8 us/call
      is_mesh_transform   2.449 ms  1.47 Mcyc   4.3% x1 2449 us/call
      mesh_visibility     1.220 ms 732.2 kcyc   2.2% x1 1220 us/call
  is_ripple_prepare       28.04 us  16.8 kcyc   0.0% x1 28.0 us/call
  canvas_clear            95.53 us  57.3 kcyc   0.2% x1 95.5 us/call
  canvas_buffer_wait       0.13 us     76 cyc   0.0% x1 0.13 us/call
```

Wall min/mean/max was 54.177/56.629/58.986 ms. Render averaged 56.629 ms per frame, including 56.324 ms in `is_timeline_step`.

This regime leaves 5.871 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 2481–2496)

```
scan_mesh_raster         33.922 ms 20.35 Mcyc  74.9% x203.44 167 us/call
scan_face_setup           3.162 ms  1.90 Mcyc   7.0% x219.69 14.4 us/call
mesh_visibility          511.77 us 307.1 kcyc   1.1% x1 512 us/call
frame                    45.294 ms 27.18 Mcyc 100.0%
  pov_preserve_half      143.39 us  86.0 kcyc   0.3% x1 143 us/call
  is_timeline_step       44.345 ms 26.61 Mcyc  97.9% x1 44345 us/call
    hk_conway_previous    1.030 ms 618.3 kcyc   2.3% x0.94 1099 us/call
    is_build_draw        41.462 ms 24.88 Mcyc  91.5% x1 41462 us/call
      is_build_scan      40.249 ms 24.15 Mcyc  88.9% x1 40249 us/call
        scan_face_setup   62.66 us  37.6 kcyc   0.1% x2.5 25.1 us/call
        mesh_copy_capture  37.59 us  22.6 kcyc   0.1% x2.62 14.3 us/call
        mesh_capture     385.67 us 231.4 kcyc   0.9% x2.81 137 us/call
        mesh_classify    354.86 us 212.9 kcyc   0.8% x0.56 631 us/call
      is_mesh_transform   62.03 us  37.2 kcyc   0.1% x1 62.0 us/call
    hk_conway_compile    411.21 us 246.7 kcyc   0.9% x1 411 us/call
    hk_conway_sweep       1.165 ms 699.3 kcyc   2.6% x1 1165 us/call
  is_ripple_prepare        7.49 us   4.5 kcyc   0.0% x1 7.49 us/call
  canvas_clear            89.17 us  53.5 kcyc   0.2% x1 89.2 us/call
  canvas_buffer_wait       0.13 us     76 cyc   0.0% x1 0.13 us/call
```

Wall min/mean/max was 41.004/45.294/50.936 ms. Render averaged 45.294 ms per frame, including 44.345 ms in `is_timeline_step`.

This regime leaves 17.206 ms of mean render headroom per display window.

### Final complete in-show window (frames 4865–4880)

```
frame                    25.678 ms 15.41 Mcyc 100.0%
  pov_preserve_half      148.27 us  89.0 kcyc   0.6% x1 148 us/call
  is_timeline_step       25.421 ms 15.25 Mcyc  99.0% x1 25421 us/call
    is_draw_shape        25.308 ms 15.18 Mcyc  98.6% x1 25308 us/call
      is_mesh_scan       24.387 ms 14.63 Mcyc  95.0% x1 24387 us/call
        scan_mesh_raster 21.125 ms 12.68 Mcyc  82.3% x173.19 122 us/call
        scan_face_setup   2.468 ms  1.48 Mcyc   9.6% x183.62 13.4 us/call
      is_face_offsets     16.55 us   9.9 kcyc   0.1% x1 16.5 us/call
      is_mesh_transform   97.34 us  58.4 kcyc   0.4% x1 97.3 us/call
      mesh_visibility    594.57 us 356.7 kcyc   2.3% x1 595 us/call
  is_ripple_prepare        6.70 us   4.0 kcyc   0.0% x1 6.70 us/call
  canvas_clear            86.18 us  51.7 kcyc   0.3% x1 86.2 us/call
  canvas_buffer_wait       0.12 us     75 cyc   0.0% x1 0.12 us/call
```

Wall min/mean/max was 25.524/25.677/26.044 ms. Render averaged 25.678 ms per frame, including 25.421 ms in `is_timeline_step`.

This regime leaves 36.822 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

Qualification requires a completed-build marker, a full-envelope window with one shape draw per frame and zero build draws. The table gives seed and built vertex/edge/face/index counts. Face visibility can change the number of raster calls within a held shape.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | completed mesh, full envelope | 9 | 2049–2064 | seed 20/30/12/60 / built 3240/4320/1082/8640 | — | 56.324 | 56.629 | 58.986 | 16 | 0/352 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | completed mesh, full envelope | 9 | 1425–1440 | seed 60/90/32/180 / built 1620/2160/542/4320 | — | 46.549 | 46.857 | 49.500 | 16 | 0/288 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | completed mesh, full envelope | 9 | 4785–4800 | seed 12/30/20/60 / built 1620/2160/542/4320 | — | 39.629 | 39.924 | 41.827 | 16 | 0/306 |
| icosahedron_snub_relax_truncate033_hankin62 | completed mesh, full envelope | 9 | 4465–4480 | seed 12/30/20/60 / built 1350/1800/452/3600 | — | 33.085 | 33.379 | 34.977 | 16 | 0/288 |
| truncatedCube_hk90_meta | completed mesh, full envelope | 9 | 1137–1152 | seed 24/36/14/72 / built 218/648/432/1296 | — | 32.430 | 32.728 | 53.457 | 16 | 0/368 |
| truncatedIcosahedron_hk58_chamfer63 | completed mesh, full envelope | 9 | 417–432 | seed 60/90/32/180 / built 990/1440/452/2880 | — | 31.494 | 31.789 | 33.263 | 16 | 0/248 |
| dodecahedron_hk72_ambo_dual_hk20 | completed mesh, full envelope | 9 | 4177–4192 | seed 20/30/12/60 / built 540/720/182/1440 | — | 30.942 | 31.249 | 34.533 | 16 | 0/352 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | completed mesh, full envelope | 9 | 2641–2656 | seed 120/180/62/360 / built 542/1080/540/2160 | — | 30.730 | 31.026 | 50.936 | 16 | 0/312 |
| rhombicuboctahedron_hk63_ambo_hk63 | completed mesh, full envelope | 9 | 3553–3568 | seed 24/48/26/96 / built 864/1152/290/2304 | — | 30.135 | 30.433 | 32.018 | 16 | 0/280 |
| icosahedron_ambo_truncate033_hankin59 | completed mesh, full envelope | 9 | 1681–1696 | seed 12/30/20/60 / built 540/720/182/1440 | — | 28.014 | 28.322 | 32.493 | 16 | 0/272 |
| snubDodecahedron_truncate5d_ambo_dual | completed mesh, full envelope | 9 | 3265–3280 | seed 60/150/92/300 / built 452/900/450/1800 | — | 26.049 | 26.346 | 34.949 | 16 | 0/312 |
| dodecahedron_hk54_ambo_hk72 | completed mesh, full envelope | 9 | 3825–3840 | seed 20/30/12/60 / built 540/720/182/1440 | — | 26.029 | 26.321 | 27.398 | 16 | 0/280 |
| dodecahedron_hk62_ambo_hk62 | completed mesh, full envelope | 9 | 177–192 | seed 20/30/12/60 / built 540/720/182/1440 | — | 25.443 | 25.738 | 26.350 | 16 | 0/278 |
| dodecahedron_bevel2_relax_gyro | completed mesh, full envelope | 9 | 769–784 | seed 20/30/12/60 / built 542/900/360/1800 | — | 25.385 | 25.680 | 34.746 | 16 | 0/352 |
| octahedron_hk17_ambo_hk73 | completed mesh, full envelope | 9 | 2321–2336 | seed 6/12/8/24 / built 216/288/74/576 | — | 25.016 | 25.327 | 36.242 | 16 | 0/280 |
| icosidodecahedron_truncate5d_ambo_dual | completed mesh, full envelope | 9 | 2929–2944 | seed 30/60/32/120 / built 182/360/180/720 | — | 19.161 | 19.463 | 27.186 | 16 | 0/312 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 2049–2064, measured over the separate 998,528 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1150.4/frame  0.465/1.776/25.613 us  3.27%
  isr_pack           143.8/frame  6.813/7.546/10.263 us  1.73%
  isr_dma_submit     143.8/frame  0.593/0.920/9.508 us  0.21%
```

Inclusive wake time leaves 60.456 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `is_timeline_step` took 56.324 ms/frame (99.5% of draw time).
2. `is_draw_shape` took 56.146 ms/frame (99.1% of draw time).
3. `is_mesh_scan` took 46.808 ms/frame (82.7% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=IslamicStars` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
