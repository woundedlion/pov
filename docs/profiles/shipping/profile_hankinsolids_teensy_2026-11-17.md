# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\hankinsolids_ship.log`. Replaces `profile_hankinsolids_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 170 s / 2720 nominal frames at 16 fps; window 16. Raw capture 210 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HankinSolids profile 173 16`; analyze only the 170 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:160448, data:162984, headers:8344; free:1699840 / RAM1: variables:314656, code:31416, padding:1352; free:176864 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 961–976 root cycles / 600 MHz = 1000056.562 us; measured wall sum 1000053 us, delta 3.56 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:22, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `72e46edfbee12186c076ab005594b88e1ebcffc06ee50746621a89a0c12ab541`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 170 s (2720 nominal frames at 16 fps). Show-frame peak render is **32.980 ms**, spilled **0/2720**. Setup frame 1 is **18.945 ms** and remains separate. Live show peak is **32.980 ms**, spilled **0/2719**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–2720. The raw log retains 3337 frame records; the 617 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 32.980 (18), spilled 🟢 0/2719 (0.00%).

README setup frame 1: 18.945 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    56.827 ms  34.096 Mcyc 100.0%
  pov_preserve_half      135.90 us   81.54 kcyc   0.2% x0.9375 145 us/call
  hk_timeline_step       17.402 ms  10.441 Mcyc  30.6%
    hk_update_hankin      25.68 us   15.41 kcyc   0.0% x0.9375 27 us/call
    hk_draw_mesh         17.131 ms  10.279 Mcyc  30.1%
      hk_mesh_scan       16.774 ms  10.065 Mcyc  29.5%
        scan_face_sweep   9.242 ms   5.545 Mcyc  16.3%
          scan_sweep_column   2.850 ms   1.710 Mcyc   5.0% x594 5 us/call
        mesh_sweep_paths   41.31 us   24.79 kcyc   0.1% x4.125 10 us/call
        scan_face_setup  147.44 us   88.46 kcyc   0.3% x4 37 us/call
        mesh_copy_capture   36.15 us   21.69 kcyc   0.1% x4 9 us/call
        mesh_capture       8.09 us    4.86 kcyc   0.0% x1.875 4 us/call
        mesh_classify     82.38 us   49.43 kcyc   0.1% x0.8125 101 us/call
        scan_mesh_raster   6.685 ms   4.011 Mcyc  11.8% x6 1114 us/call
        scan_face_setup   63.99 us   38.39 kcyc   0.1% x2.25 28 us/call
      hk_mesh_transform  160.80 us   96.48 kcyc   0.3% x1 161 us/call
    hk_update_hankin      76.13 us   45.68 kcyc   0.1% x1 76 us/call
  canvas_clear            98.26 us   58.95 kcyc   0.2% x1 98 us/call
  canvas_buffer_wait     39.135 ms  23.481 Mcyc  68.9% x1 39135 us/call
```

Wall min/avg/max = 16.536/56.827/64.892 ms. `hk_timeline_step` averages 17.402 ms/f; total render averages 17.692 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 961–976)

```
frame                    62.504 ms  37.502 Mcyc 100.0%
  pov_preserve_half      144.72 us   86.83 kcyc   0.2% x1 145 us/call
  hk_timeline_step       28.880 ms  17.328 Mcyc  46.2%
    hk_update_hankin     132.11 us   79.27 kcyc   0.2% x1 132 us/call
    hk_draw_mesh         28.358 ms  17.015 Mcyc  45.4%
      hk_mesh_scan       27.830 ms  16.698 Mcyc  44.5%
        scan_face_sweep   7.520 ms   4.512 Mcyc  12.0%
          scan_sweep_column   2.102 ms   1.261 Mcyc   3.4% x290.75 7 us/call
        mesh_sweep_paths   33.75 us   20.25 kcyc   0.1% x2.25 15 us/call
        scan_face_setup  168.98 us  101.39 kcyc   0.3% x10.6875 16 us/call
        mesh_copy_capture   37.34 us   22.40 kcyc   0.1% x11.125 3 us/call
        mesh_capture      86.88 us   52.13 kcyc   0.1% x2.875 30 us/call
        mesh_classify    145.03 us   87.02 kcyc   0.2% x0.875 166 us/call
        scan_mesh_raster  18.486 ms  11.092 Mcyc  29.6% x27.1875 680 us/call
        scan_face_setup  472.74 us  283.64 kcyc   0.8% x21 23 us/call
      hk_mesh_transform  311.46 us  186.88 kcyc   0.5% x1 311 us/call
    hk_update_hankin     204.51 us  122.70 kcyc   0.3% x1 205 us/call
  canvas_clear            89.51 us   53.71 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     33.314 ms  19.988 Mcyc  53.3% x1 33314 us/call
```

Wall min/avg/max = 53.441/62.503/69.364 ms. `hk_timeline_step` averages 28.880 ms/f; total render averages 29.189 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 1409–1424)

```
frame                    62.276 ms  37.366 Mcyc 100.0%
  pov_preserve_half      146.90 us   88.14 kcyc   0.2% x1 147 us/call
  hk_timeline_step       29.813 ms  17.888 Mcyc  47.9%
    hk_update_hankin     257.01 us  154.20 kcyc   0.4% x1 257 us/call
    hk_draw_mesh         29.059 ms  17.436 Mcyc  46.7%
      hk_mesh_scan       28.420 ms  17.052 Mcyc  45.6%
        scan_face_sweep   2.308 ms   1.385 Mcyc   3.7%
          scan_sweep_column  747.71 us  448.63 kcyc   1.2% x72 10 us/call
        mesh_sweep_paths    8.23 us    4.94 kcyc   0.0% x0.5 16 us/call
        scan_face_setup   49.03 us   29.42 kcyc   0.1% x1.875 26 us/call
        mesh_copy_capture   19.37 us   11.62 kcyc   0.0% x1.875 10 us/call
        mesh_capture     118.99 us   71.39 kcyc   0.2% x3.75 32 us/call
        mesh_classify    131.35 us   78.81 kcyc   0.2% x0.75 175 us/call
        scan_mesh_raster  24.311 ms  14.587 Mcyc  39.0% x44.75 543 us/call
        scan_face_setup  757.12 us  454.27 kcyc   1.2% x46.375 16 us/call
      hk_mesh_transform  393.65 us  236.19 kcyc   0.6% x1 394 us/call
    hk_update_hankin     324.37 us  194.62 kcyc   0.5% x1 324 us/call
  canvas_clear            89.28 us   53.57 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     32.164 ms  19.299 Mcyc  51.6% x1 32164 us/call
```

Wall min/avg/max = 55.998/62.276/68.728 ms. `hk_timeline_step` averages 29.813 ms/f; total render averages 30.112 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 2161–2176)

```
frame                    62.930 ms  37.758 Mcyc 100.0%
  pov_preserve_half      143.04 us   85.82 kcyc   0.2% x1 143 us/call
  hk_timeline_step       20.895 ms  12.537 Mcyc  33.2%
    hk_conway_previous    22.36 us   13.41 kcyc   0.0% x0.1875 119 us/call
    hk_draw_mesh          2.105 ms   1.263 Mcyc   3.3%
      hk_mesh_scan        2.071 ms   1.243 Mcyc   3.3% x0.1875 11048 us/call
      hk_mesh_transform    1.09 us    0.65 kcyc   0.0% x0.1875 6 us/call
    hk_conway_compile      9.00 us    5.40 kcyc   0.0% x0.1875 48 us/call
    hk_conway_sweep       31.27 us   18.76 kcyc   0.0% x0.1875 167 us/call
    hk_update_hankin     294.17 us  176.50 kcyc   0.5% x0.75 392 us/call
    hk_draw_mesh         17.690 ms  10.614 Mcyc  28.1%
      hk_mesh_scan       17.452 ms  10.471 Mcyc  27.7%
        scan_mesh_raster  18.373 ms  11.024 Mcyc  29.2% x46.625 394 us/call
        scan_face_setup  663.79 us  398.27 kcyc   1.1% x50.0625 13 us/call
      hk_mesh_transform   23.65 us   14.19 kcyc   0.0% x0.8125 29 us/call
    hk_update_hankin     360.97 us  216.58 kcyc   0.6% x0.8125 444 us/call
  canvas_clear            88.76 us   53.26 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     41.760 ms  25.056 Mcyc  66.4% x1 41760 us/call
```

Wall min/avg/max = 60.222/62.929/70.606 ms. `hk_timeline_step` averages 20.895 ms/f; total render averages 21.169 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 2705–2720)

```
frame                    62.519 ms  37.511 Mcyc 100.0%
  pov_preserve_half      142.82 us   85.69 kcyc   0.2% x1 143 us/call
  hk_timeline_step       16.594 ms   9.956 Mcyc  26.5%
    hk_conway_previous     3.58 us    2.15 kcyc   0.0% x0.125 29 us/call
    hk_draw_mesh          2.049 ms   1.230 Mcyc   3.3%
      hk_mesh_scan        2.026 ms   1.216 Mcyc   3.2% x0.125 16210 us/call
      hk_mesh_transform    0.92 us    0.55 kcyc   0.0% x0.125 7 us/call
    hk_conway_compile      4.01 us    2.41 kcyc   0.0% x0.125 32 us/call
    hk_conway_sweep        9.47 us    5.68 kcyc   0.0% x0.125 76 us/call
    hk_update_hankin      21.68 us   13.01 kcyc   0.0% x0.8125 27 us/call
    hk_draw_mesh         14.223 ms   8.534 Mcyc  22.8%
      hk_mesh_scan       13.916 ms   8.349 Mcyc  22.3%
        scan_face_sweep   6.032 ms   3.619 Mcyc   9.6%
          scan_sweep_column   1.704 ms   1.022 Mcyc   2.7% x342 5 us/call
        mesh_sweep_paths   23.14 us   13.89 kcyc   0.0% x2.375 10 us/call
        scan_face_setup  114.66 us   68.80 kcyc   0.2% x3.875 30 us/call
        mesh_copy_capture   33.39 us   20.03 kcyc   0.1% x3.875 9 us/call
        mesh_capture      16.26 us    9.75 kcyc   0.0% x1.875 9 us/call
        mesh_classify     70.26 us   42.16 kcyc   0.1% x0.6875 102 us/call
        scan_mesh_raster   8.987 ms   5.392 Mcyc  14.4% x8.625 1042 us/call
        scan_face_setup  148.67 us   89.20 kcyc   0.2% x5.5 27 us/call
      hk_mesh_transform  136.09 us   81.66 kcyc   0.2% x0.875 156 us/call
    hk_update_hankin      68.56 us   41.14 kcyc   0.1% x0.875 78 us/call
  canvas_clear            89.77 us   53.86 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     45.622 ms  27.373 Mcyc  73.0% x1 45622 us/call
```

Wall min/avg/max = 56.334/62.518/68.648 ms. `hk_timeline_step` averages 16.594 ms/f; total render averages 16.896 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Qualified windows contain one solid and one hk_update_hankin call per frame. This selects complete fixed-solid Hankin-angle sweeps and excludes Conway construction legs. Face raster counts vary with visibility, so no modal raster-call filter is applied. The angle, palette and strap handoffs remain live.

The show-scoped data covers 18 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| truncatedCuboctahedron | 4 | 1409–1424 | fixed-solid Hankin-angle sweep | — | — | 29.813 | 30.112 | 16 | 32.590 | 0/117 |
| truncatedCube | 5 | 961–976 | fixed-solid Hankin-angle sweep | — | — | 28.880 | 29.189 | 16 | 32.980 | 0/105 |
| snubDodecahedron | 4 | 2545–2560 | fixed-solid Hankin-angle sweep | — | — | 27.904 | 28.191 | 16 | 32.091 | 0/117 |
| rhombicosidodecahedron | 4 | 2193–2208 | fixed-solid Hankin-angle sweep | — | — | 27.214 | 27.499 | 16 | 30.456 | 0/117 |
| truncatedIcosahedron | 4 | 2033–2048 | fixed-solid Hankin-angle sweep | — | — | 26.747 | 27.028 | 16 | 28.997 | 0/105 |
| truncatedIcosidodecahedron | 4 | 2369–2384 | fixed-solid Hankin-angle sweep | — | — | 26.650 | 26.946 | 16 | 30.452 | 0/117 |
| truncatedDodecahedron | 4 | 1905–1920 | fixed-solid Hankin-angle sweep | — | — | 26.593 | 26.896 | 16 | 32.318 | 0/105 |
| truncatedOctahedron | 5 | 1089–1104 | fixed-solid Hankin-angle sweep | — | — | 24.327 | 24.631 | 16 | 26.520 | 0/105 |
| snubCube | 5 | 1585–1600 | fixed-solid Hankin-angle sweep | — | — | 23.429 | 23.714 | 16 | 29.619 | 0/117 |
| truncatedTetrahedron | 4 | 737–752 | fixed-solid Hankin-angle sweep | — | — | 23.137 | 23.451 | 16 | 28.641 | 0/105 |
| icosidodecahedron | 4 | 1793–1808 | fixed-solid Hankin-angle sweep | — | — | 22.937 | 23.216 | 16 | 27.178 | 0/285 |
| rhombicuboctahedron | 4 | 1249–1264 | fixed-solid Hankin-angle sweep | — | — | 20.996 | 21.281 | 16 | 24.074 | 0/105 |
| dodecahedron | 4 | 465–480 | fixed-solid Hankin-angle sweep | — | — | 20.105 | 20.393 | 16 | 23.120 | 0/225 |
| octahedron | 4 | 289–304 | fixed-solid Hankin-angle sweep | — | — | 19.932 | 20.239 | 16 | 22.621 | 0/177 |
| cuboctahedron | 4 | 833–848 | fixed-solid Hankin-angle sweep | — | — | 19.120 | 19.427 | 16 | 21.821 | 0/285 |
| cube | 4 | 193–208 | fixed-solid Hankin-angle sweep | — | — | 18.845 | 19.143 | 16 | 23.662 | 0/213 |
| icosahedron | 4 | 625–640 | fixed-solid Hankin-angle sweep | — | — | 18.602 | 18.926 | 16 | 23.192 | 0/177 |
| tetrahedron | 4 | 17–32 | fixed-solid Hankin-angle sweep | — | — | 16.145 | 16.469 | 16 | 19.323 | 0/142 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 961–976; CPU shares use its separately snapshotted ISR interval of 1002554 us.

```
isr_wake       1155.0/f 0.288/1.703/29.728 us 3.14%
  isr_pack       144.4/f 7.780/8.591/10.990 us 1.97%
  isr_dma_submit 144.4/f 0.680/0.918/7.651 us 0.21%
```

- Mean pack 8.591 us/call versus submit 0.918 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.14%, corresponding to 60.538 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.283×; Peak live show render: 0.467×; Costliest qualified preset: 0.482×; Transition / non-hold regime: 0.339×; Show end: 0.270×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `hk_timeline_step` — 46.2% of the frame, 28.880 ms/f.
2. `hk_draw_mesh` — 45.4% of the frame, 28.358 ms/f.
3. `hk_mesh_scan` — 44.5% of the frame, 27.830 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HankinSolids`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
