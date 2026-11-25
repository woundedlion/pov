# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 170-second Phantasm slot. Raw capture `C:/work/temp/holosphere-show-sweep-6cca5fb7/hankinsolids_ship.log`. Supersedes `profile_hankinsolids_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | HankinSolids at source `6cca5fb7f88cbdfb7022d777aa93be143fb957be` |
| Method | 170 s slot, 176 s acquisition, 16 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HankinSolids profile 176 16` from the captured revision |

Captured 2026-11-20 15:43, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `9a16ffe5627466cbb5e89888da79b8f643f75f1de87fcaf85a3d53b9392c8c70`.

Profile image size `FLASH: code:160976, data:163008, headers:8816 free for files:1698816 / RAM1: variables:314656, code:31256, padding:1512 free for local variables:176864 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111904, data:724920, headers:8400 free for files:186392 / RAM1: variables:313760, code:164856, padding:31752 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `7e91d369ccfba4019b40d947e6511528ad0aaa6d0e68e4e9abb4f5cb04481928`. Shipping Phantasm ELF SHA-256 `a3a44841c360a2ca758ee00eae95793639ed0d7d7c26b1fc4944653f463d7f95`.

Exactness cross-check, frames 1409–1424. The root counter measured 481932.658 µs, compared with 481,935 µs of draw time (4.86 ppm).

## Frame cadence

**Peak confirmed live render 32.517 ms** at frame 1420 (truncatedCuboctahedron). **0/2,703 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 19.007 ms. Live render averaged 17.310 ms.

Among complete live readouts, the largest window mean was 30.119 ms/frame at frames 1409–1424. The dominant scope averaged 17.037 ms/frame across those readouts.

The shipping ELF assigns 1,360 revolutions (170 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 170 s. Frames 2–2704 are confirmed live. The raw file contains 2,793 frames, including 73 confirmed after the boundary.

Boundary resolution leaves 16 frame(s), 2705–2720, with uncertain membership. Their largest render was 20.268 ms and 0 exceeded the budget. The live denominator lies between 2,703 and 2,719. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 32.517 (18), spilled 🟢 0/2703 (0.00%).

README setup frame 1: 19.007 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–16)

```
frame                    17.233 ms 10.34 Mcyc 100.0%
  pov_preserve_half      134.75 us  80.8 kcyc   0.8% x0.94 144 us/call
  hk_timeline_step       16.967 ms 10.18 Mcyc  98.5% x1 16967 us/call
    hk_update_hankin      23.72 us  14.2 kcyc   0.1% x0.94 25.3 us/call
    hk_draw_mesh         16.729 ms 10.04 Mcyc  97.1% x1 16729 us/call
      hk_mesh_scan       16.363 ms  9.82 Mcyc  95.0% x1 16363 us/call
        scan_face_sweep   8.849 ms  5.31 Mcyc  51.3% x4.12 2145 us/call
          scan_sweep_column  2.632 ms  1.58 Mcyc  15.3% x594 4.43 us/call
        mesh_sweep_paths  37.95 us  22.8 kcyc   0.2% x4.12 9.20 us/call
        scan_face_setup  144.82 us  86.9 kcyc   0.8% x4 36.2 us/call
        mesh_copy_capture  44.21 us  26.5 kcyc   0.3% x4 11.1 us/call
        mesh_capture       8.73 us   5.2 kcyc   0.1% x1.88 4.65 us/call
        mesh_classify     80.78 us  48.5 kcyc   0.5% x0.81 99.4 us/call
        scan_mesh_raster  6.681 ms  4.01 Mcyc  38.8% x6 1113 us/call
        scan_face_setup   67.48 us  40.5 kcyc   0.4% x2.25 30.0 us/call
      hk_mesh_transform  141.16 us  84.7 kcyc   0.8% x1 141 us/call
    hk_update_hankin      71.17 us  42.7 kcyc   0.4% x1 71.2 us/call
  canvas_clear           100.69 us  60.4 kcyc   0.6% x1 101 us/call
  canvas_buffer_wait       1.40 us    837 cyc   0.0% x1 1.40 us/call
```

Wall min/mean/max was 15.555/17.233/19.007 ms. Render averaged 17.232 ms per frame, including 16.967 ms in `hk_timeline_step`.

This regime leaves 45.269 ms of mean render headroom per display window.

### Window containing the live peak (frames 1409–1424)

```
frame                    30.121 ms 18.07 Mcyc 100.0%
  pov_preserve_half      144.89 us  86.9 kcyc   0.5% x1 145 us/call
  hk_timeline_step       29.848 ms 17.91 Mcyc  99.1% x1 29848 us/call
    hk_update_hankin     255.14 us 153.1 kcyc   0.8% x1 255 us/call
    hk_draw_mesh         29.128 ms 17.48 Mcyc  96.7% x1 29128 us/call
      hk_mesh_scan       28.474 ms 17.08 Mcyc  94.5% x1 28474 us/call
        scan_face_sweep   2.241 ms  1.34 Mcyc   7.4% x0.5 4481 us/call
          scan_sweep_column 720.40 us 432.2 kcyc   2.4% x72 10.0 us/call
        mesh_sweep_paths   7.62 us   4.6 kcyc   0.0% x0.5 15.2 us/call
        scan_face_setup   43.76 us  26.3 kcyc   0.1% x1.88 23.3 us/call
        mesh_copy_capture  15.40 us   9.2 kcyc   0.1% x1.88 8.22 us/call
        mesh_capture     127.15 us  76.3 kcyc   0.4% x3.75 33.9 us/call
        mesh_classify    125.67 us  75.4 kcyc   0.4% x0.75 168 us/call
        scan_mesh_raster 24.471 ms 14.68 Mcyc  81.2% x44.75 547 us/call
        scan_face_setup  767.04 us 460.2 kcyc   2.5% x46.38 16.5 us/call
      hk_mesh_transform  378.23 us 226.9 kcyc   1.3% x1 378 us/call
    hk_update_hankin     311.94 us 187.2 kcyc   1.0% x1 312 us/call
  canvas_clear            91.17 us  54.7 kcyc   0.3% x1 91.2 us/call
  canvas_buffer_wait       1.33 us    796 cyc   0.0% x1 1.33 us/call
```

Wall min/mean/max was 23.960/30.120/32.518 ms. Render averaged 30.119 ms per frame, including 29.848 ms in `hk_timeline_step`.

This regime leaves 32.381 ms of mean render headroom per display window.

### Highest render outside qualified holds (frames 2161–2176)

```
frame                    21.071 ms 12.64 Mcyc 100.0%
  pov_preserve_half      144.41 us  86.6 kcyc   0.7% x1 144 us/call
  hk_timeline_step       20.792 ms 12.48 Mcyc  98.7% x1 20792 us/call
    hk_conway_previous    19.91 us  11.9 kcyc   0.1% x0.19 106 us/call
    hk_draw_mesh          2.091 ms  1.25 Mcyc   9.9% x0.19 11150 us/call
      hk_mesh_scan        2.057 ms  1.23 Mcyc   9.8% x0.19 10972 us/call
      hk_mesh_transform    0.85 us    509 cyc   0.0% x0.19 4.52 us/call
    hk_conway_compile      8.57 us   5.1 kcyc   0.0% x0.19 45.7 us/call
    hk_conway_sweep       31.53 us  18.9 kcyc   0.1% x0.19 168 us/call
    hk_update_hankin     294.46 us 176.7 kcyc   1.4% x0.75 393 us/call
    hk_draw_mesh         17.640 ms 10.58 Mcyc  83.7% x0.81 21710 us/call
      hk_mesh_scan       17.381 ms 10.43 Mcyc  82.5% x0.81 21392 us/call
        scan_mesh_raster 18.336 ms 11.00 Mcyc  87.0% x46.62 393 us/call
        scan_face_setup  648.01 us 388.8 kcyc   3.1% x50.06 12.9 us/call
      hk_mesh_transform   17.93 us  10.8 kcyc   0.1% x0.81 22.1 us/call
    hk_update_hankin     348.88 us 209.3 kcyc   1.7% x0.81 429 us/call
  canvas_clear            90.69 us  54.4 kcyc   0.4% x1 90.7 us/call
  canvas_buffer_wait       1.36 us    814 cyc   0.0% x1 1.36 us/call
```

Wall min/mean/max was 11.700/21.070/26.643 ms. Render averaged 21.069 ms per frame, including 20.792 ms in `hk_timeline_step`.

This regime leaves 41.431 ms of mean render headroom per display window.

### Final complete in-show window (frames 2689–2704)

```
scan_face_sweep           4.654 ms  2.79 Mcyc  26.2% x2.25 2068 us/call
  scan_sweep_column       1.125 ms 675.0 kcyc   6.3% x306.56 3.67 us/call
mesh_sweep_paths          17.06 us  10.2 kcyc   0.1% x2.25 7.58 us/call
scan_face_setup           87.49 us  52.5 kcyc   0.5% x6.31 13.9 us/call
mesh_copy_capture         19.35 us  11.6 kcyc   0.1% x6.38 3.03 us/call
mesh_capture              80.76 us  48.5 kcyc   0.5% x3.19 25.3 us/call
mesh_classify             98.32 us  59.0 kcyc   0.6% x0.94 105 us/call
scan_mesh_raster         10.711 ms  6.43 Mcyc  60.3% x18.25 587 us/call
scan_face_setup          327.54 us 196.5 kcyc   1.8% x14.44 22.7 us/call
frame                    17.761 ms 10.66 Mcyc 100.0%
  pov_preserve_half      146.35 us  87.8 kcyc   0.8% x1 146 us/call
  hk_timeline_step       17.488 ms 10.49 Mcyc  98.5% x1 17488 us/call
    hk_conway_previous    34.24 us  20.5 kcyc   0.2% x1 34.2 us/call
    hk_draw_mesh         17.088 ms 10.25 Mcyc  96.2% x1 17088 us/call
      hk_mesh_scan       16.662 ms 10.00 Mcyc  93.8% x1 16662 us/call
      hk_mesh_transform  211.74 us 127.0 kcyc   1.2% x1 212 us/call
    hk_conway_compile     46.57 us  27.9 kcyc   0.3% x1 46.6 us/call
    hk_conway_sweep      125.21 us  75.1 kcyc   0.7% x1 125 us/call
  canvas_clear            92.14 us  55.3 kcyc   0.5% x1 92.1 us/call
  canvas_buffer_wait       1.31 us    784 cyc   0.0% x1 1.31 us/call
```

Wall min/mean/max was 13.585/17.761/22.187 ms. Render averaged 17.760 ms per frame, including 17.488 ms in `hk_timeline_step`.

This regime leaves 44.740 ms of mean render headroom per display window.

### Per-preset table

Live peaks and spills include transitions attributed to each serial owner. Each row selects the eligible complete window with the largest dominant-scope cost. Bounds exclude setup and mixed-owner windows. The fps column gives the quantized cadence implied by the worst live render for that owner.

| Preset / shape | Regime | Eligible windows | Frames | V/E/F/I | Blends/frame | Scope ms/frame | Render ms/frame | Live peak ms | fps | Live spilled |
|---|---|--:|---|---|--:|--:|--:|--:|--:|--:|
| truncatedCuboctahedron | fixed-solid angle sweep | 4 | 1409–1424 | — | — | 29.848 | 30.119 | 32.517 | 16 | 0/117 |
| truncatedCube | fixed-solid angle sweep | 5 | 961–976 | — | — | 28.564 | 28.834 | 32.515 | 16 | 0/105 |
| snubDodecahedron | fixed-solid angle sweep | 4 | 2545–2560 | — | — | 27.917 | 28.187 | 32.118 | 16 | 0/117 |
| rhombicosidodecahedron | fixed-solid angle sweep | 4 | 2193–2208 | — | — | 27.193 | 27.465 | 30.394 | 16 | 0/117 |
| truncatedIcosahedron | fixed-solid angle sweep | 4 | 2033–2048 | — | — | 26.900 | 27.170 | 29.146 | 16 | 0/105 |
| truncatedIcosidodecahedron | fixed-solid angle sweep | 4 | 2369–2384 | — | — | 26.773 | 27.040 | 30.618 | 16 | 0/117 |
| truncatedDodecahedron | fixed-solid angle sweep | 4 | 1905–1920 | — | — | 26.400 | 26.668 | 31.690 | 16 | 0/105 |
| truncatedOctahedron | fixed-solid angle sweep | 5 | 1089–1104 | — | — | 24.321 | 24.588 | 26.695 | 16 | 0/105 |
| snubCube | fixed-solid angle sweep | 5 | 1553–1568 | — | — | 23.369 | 23.635 | 29.497 | 16 | 0/117 |
| icosidodecahedron | fixed-solid angle sweep | 4 | 1793–1808 | — | — | 22.978 | 23.249 | 27.246 | 16 | 0/285 |
| truncatedTetrahedron | fixed-solid angle sweep | 4 | 737–752 | — | — | 22.548 | 22.817 | 28.362 | 16 | 0/105 |
| rhombicuboctahedron | fixed-solid angle sweep | 4 | 1249–1264 | — | — | 20.918 | 21.193 | 23.984 | 16 | 0/105 |
| dodecahedron | fixed-solid angle sweep | 4 | 465–480 | — | — | 20.193 | 20.462 | 23.292 | 16 | 0/225 |
| octahedron | fixed-solid angle sweep | 4 | 289–304 | — | — | 19.285 | 19.553 | 22.212 | 16 | 0/177 |
| cuboctahedron | fixed-solid angle sweep | 4 | 849–864 | — | — | 18.867 | 19.136 | 21.639 | 16 | 0/285 |
| cube | fixed-solid angle sweep | 4 | 193–208 | — | — | 18.578 | 18.851 | 23.223 | 16 | 0/213 |
| icosahedron | fixed-solid angle sweep | 4 | 625–640 | — | — | 18.528 | 18.797 | 23.335 | 16 | 0/176 |
| tetrahedron | fixed-solid angle sweep | 4 | 17–32 | — | — | 15.663 | 15.934 | 18.406 | 16 | 0/127 |

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 1409–1424, measured over the separate 998,731 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1150.6/frame  0.498/1.720/13.873 us  3.17%
  isr_pack           143.8/frame  6.760/7.298/9.943 us  1.68%
  isr_dma_submit     143.8/frame  0.585/0.915/1.673 us  0.21%
```

Inclusive wake time leaves 60.519 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `hk_timeline_step` took 29.848 ms/frame (99.1% of draw time).
2. `hk_draw_mesh` took 29.128 ms/frame (96.7% of draw time).
3. `hk_mesh_scan` took 28.474 ms/frame (94.5% of draw time).

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. Scoped optimization regions include `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=HankinSolids` and `HS_PROFILE_WINDOW=16`. Both the measured checkout and recorded source diff were clean.
