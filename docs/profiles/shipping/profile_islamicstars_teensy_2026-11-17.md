# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\islamicstars_ship.log`. Replaces `profile_islamicstars_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 306 s / 4896 nominal frames at 16 fps; window 16. Raw capture 321 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh IslamicStars profile 309 16`; analyze only the 306 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:185400, data:196164, headers:8580; free:1641472 / RAM1: variables:314688, code:30680, padding:2088; free:176832 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 2049–2064 root cycles / 600 MHz = 995316.953 us; measured wall sum 995317 us, delta 0.05 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:09, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `1f3f3c29dd08bb62ea71877f8851da7b5108ebc93403a21fd94ee34a212b3d47`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 306 s (4896 nominal frames at 16 fps). Show-frame peak render is **58.059 ms**, spilled **0/4896**. Setup frame 1 is **16.537 ms** and remains separate. Live show peak is **58.059 ms**, spilled **0/4895**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–4896. The raw log retains 5112 frame records; the 216 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 58.059 (16), spilled 🟢 0/4895 (0.00%).

README setup frame 1: 16.537 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    56.314 ms  33.788 Mcyc 100.0%
  pov_preserve_half      134.41 us   80.64 kcyc   0.2% x0.9375 143 us/call
  is_timeline_step       10.166 ms   6.099 Mcyc  18.1%
    is_draw_shape        10.033 ms   6.020 Mcyc  17.8%
      is_mesh_scan        9.761 ms   5.857 Mcyc  17.3%
        scan_mesh_raster   9.411 ms   5.646 Mcyc  16.7% x7.9375 1186 us/call
        scan_face_setup  138.06 us   82.84 kcyc   0.2% x7.9375 17 us/call
      is_face_offsets     21.65 us   12.99 kcyc   0.0% x1 22 us/call
      is_mesh_transform   23.46 us   14.08 kcyc   0.0% x1 23 us/call
      mesh_visibility     23.24 us   13.95 kcyc   0.0% x1 23 us/call
  is_ripple_prepare        5.57 us    3.34 kcyc   0.0% x1 6 us/call
  canvas_clear            96.83 us   58.10 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     45.869 ms  27.522 Mcyc  81.5% x1 45869 us/call
```

Wall min/avg/max = 9.342/56.314/63.751 ms. `is_timeline_step` averages 10.166 ms/f; total render averages 10.445 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 2049–2064)

```
scan_face_setup          345.25 us  207.15 kcyc   0.6% x17.875 19 us/call
mesh_copy_capture        236.57 us  141.94 kcyc   0.4% x22.1875 11 us/call
mesh_capture              2.666 ms   1.600 Mcyc   4.3% x5 533 us/call
mesh_classify             1.502 ms  901.19 kcyc   2.4% x1 1502 us/call
frame                    62.207 ms  37.324 Mcyc 100.0%
  pov_preserve_half      141.66 us   85.00 kcyc   0.2% x1 142 us/call
  is_timeline_step       55.539 ms  33.323 Mcyc  89.3%
    is_draw_shape        55.384 ms  33.231 Mcyc  89.0%
      is_mesh_scan       46.091 ms  27.655 Mcyc  74.1%
        scan_face_sweep   5.529 ms   3.318 Mcyc   8.9%
          scan_sweep_column 2.342 ms 1.405 Mcyc 3.8% x143.25 16 us/call
        mesh_sweep_paths  115.71 us   69.43 kcyc   0.2% x1 116 us/call
        scan_mesh_raster  26.916 ms  16.149 Mcyc  43.3% x387.875 69 us/call
        scan_face_setup   5.501 ms   3.301 Mcyc   8.8% x426.875 13 us/call
      is_face_offsets     45.11 us   27.07 kcyc   0.1% x1 45 us/call
      is_mesh_transform   2.459 ms   1.476 Mcyc   4.0% x1 2459 us/call
      mesh_visibility     1.228 ms  736.54 kcyc   2.0% x1 1228 us/call
  is_ripple_prepare       21.33 us   12.80 kcyc   0.0% x1 21 us/call
  canvas_clear            93.54 us   56.12 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait      6.360 ms   3.816 Mcyc  10.2% x1 6360 us/call
```

Wall min/avg/max = 59.186/62.207/65.270 ms. `is_timeline_step` averages 55.539 ms/f; total render averages 55.847 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 2481–2496)

```
scan_mesh_raster         34.064 ms  20.439 Mcyc  54.3% x203.438 167 us/call
scan_face_setup           3.087 ms   1.852 Mcyc   4.9% x219.688 14 us/call
mesh_visibility          526.62 us  315.97 kcyc   0.8% x1 527 us/call
frame                    62.787 ms  37.672 Mcyc 100.0%
  pov_preserve_half      139.52 us   83.71 kcyc   0.2% x1 140 us/call
  is_timeline_step       44.151 ms  26.491 Mcyc  70.3%
    hk_conway_previous    1.032 ms  619.12 kcyc   1.6% x0.9375 1101 us/call
    is_build_draw        41.370 ms  24.822 Mcyc  65.9%
      is_build_scan      40.127 ms  24.076 Mcyc  63.9%
        scan_face_setup   49.82 us   29.89 kcyc   0.1% x2.5 20 us/call
        mesh_copy_capture   21.37 us   12.82 kcyc   0.0% x2.625 8 us/call
        mesh_capture     354.25 us  212.55 kcyc   0.6% x2.8125 126 us/call
        mesh_classify    335.15 us  201.09 kcyc   0.5% x0.5625 596 us/call
      is_mesh_transform   84.66 us   50.80 kcyc   0.1% x1 85 us/call
    hk_conway_compile    389.14 us  233.48 kcyc   0.6% x1 389 us/call
    hk_conway_sweep       1.154 ms  692.22 kcyc   1.8% x1 1154 us/call
  is_ripple_prepare        6.70 us    4.02 kcyc   0.0% x1 7 us/call
  canvas_clear            90.64 us   54.38 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     17.684 ms  10.610 Mcyc  28.2% x1 17684 us/call
```

Wall min/avg/max = 52.615/62.787/74.341 ms. `is_timeline_step` averages 44.151 ms/f; total render averages 45.104 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 4881–4896)

```
frame                    62.476 ms  37.485 Mcyc 100.0%
  pov_preserve_half      148.88 us   89.33 kcyc   0.2% x1 149 us/call
  is_timeline_step       26.042 ms  15.625 Mcyc  41.7%
    is_draw_shape        25.873 ms  15.524 Mcyc  41.4%
      is_mesh_scan       24.381 ms  14.628 Mcyc  39.0%
        scan_mesh_raster  21.018 ms  12.611 Mcyc  33.6% x180.625 116 us/call
        scan_face_setup   2.534 ms   1.520 Mcyc   4.1% x191 13 us/call
      is_face_offsets    305.68 us  183.41 kcyc   0.5% x1 306 us/call
      is_mesh_transform  101.44 us   60.87 kcyc   0.2% x1 101 us/call
      mesh_visibility    604.61 us  362.77 kcyc   1.0% x1 605 us/call
  is_ripple_prepare        4.86 us    2.92 kcyc   0.0% x1 5 us/call
  canvas_clear            88.79 us   53.28 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     36.148 ms  21.689 Mcyc  57.9% x1 36148 us/call
```

Wall min/avg/max = 61.360/62.475/63.685 ms. `is_timeline_step` averages 26.042 ms/f; total render averages 26.328 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Qualified windows have one shape owner, one is_draw_shape call per frame, no is_build_draw calls, and lie after the source-derived build/fade-in interval and before fade-out. The current visibility pass culls faces before scan_mesh_raster, so that counter varies with orientation and is not equal to the complete face count. Seed and completed geometry are listed separately. Completed geometry remains animated by the shipping ripples and orientation.

The show-scoped data covers 16 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 9 | 2049–2064 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 3240/4320/1082/8640 | — | 55.539 | 55.847 | 16 | 58.059 | 0/352 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 9 | 1425–1440 | built-shape full-envelope live-ripple window | seed 60/90/32/180; built 1620/2160/542/4320 | — | 45.912 | 46.215 | 16 | 48.857 | 0/288 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 9 | 4785–4800 | built-shape full-envelope live-ripple window | seed 12/30/20/60; built 1620/2160/542/4320 | — | 38.901 | 39.204 | 16 | 41.159 | 0/320 |
| icosahedron_snub_relax_truncate033_hankin62 | 9 | 4465–4480 | built-shape full-envelope live-ripple window | seed 12/30/20/60; built 1350/1800/452/3600 | — | 32.411 | 32.714 | 16 | 34.276 | 0/288 |
| truncatedCube_hk90_meta | 9 | 1137–1152 | built-shape full-envelope live-ripple window | seed 24/36/14/72; built 218/648/432/1296 | — | 32.061 | 32.367 | 16 | 52.982 | 0/368 |
| truncatedIcosahedron_hk58_chamfer63 | 9 | 417–432 | built-shape full-envelope live-ripple window | seed 60/90/32/180; built 990/1440/452/2880 | — | 31.079 | 31.381 | 16 | 32.937 | 0/248 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 9 | 2641–2656 | built-shape full-envelope live-ripple window | seed 120/180/62/360; built 542/1080/540/2160 | — | 30.590 | 30.896 | 16 | 51.005 | 0/312 |
| dodecahedron_hk72_ambo_dual_hk20 | 9 | 4177–4192 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 540/720/182/1440 | — | 30.471 | 30.784 | 16 | 34.073 | 0/352 |
| rhombicuboctahedron_hk63_ambo_hk63 | 9 | 3553–3568 | built-shape full-envelope live-ripple window | seed 24/48/26/96; built 864/1152/290/2304 | — | 29.476 | 29.785 | 16 | 31.533 | 0/280 |
| icosahedron_ambo_truncate033_hankin59 | 9 | 1681–1696 | built-shape full-envelope live-ripple window | seed 12/30/20/60; built 540/720/182/1440 | — | 27.536 | 27.843 | 16 | 31.911 | 0/272 |
| snubDodecahedron_truncate5d_ambo_dual | 9 | 3265–3280 | built-shape full-envelope live-ripple window | seed 60/150/92/300; built 452/900/450/1800 | — | 25.702 | 26.008 | 16 | 34.770 | 0/312 |
| dodecahedron_hk54_ambo_hk72 | 9 | 3825–3840 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 540/720/182/1440 | — | 25.482 | 25.790 | 16 | 26.867 | 0/280 |
| dodecahedron_bevel2_relax_gyro | 9 | 769–784 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 542/900/360/1800 | — | 25.109 | 25.413 | 16 | 34.551 | 0/352 |
| dodecahedron_hk62_ambo_hk62 | 9 | 177–192 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 540/720/182/1440 | — | 24.804 | 25.112 | 16 | 26.714 | 0/279 |
| octahedron_hk17_ambo_hk73 | 9 | 2321–2336 | built-shape full-envelope live-ripple window | seed 6/12/8/24; built 216/288/74/576 | — | 24.574 | 24.883 | 16 | 35.737 | 0/280 |
| icosidodecahedron_truncate5d_ambo_dual | 9 | 2929–2944 | built-shape full-envelope live-ripple window | seed 30/60/32/120; built 182/360/180/720 | — | 18.980 | 19.280 | 16 | 27.005 | 0/312 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 2049–2064; CPU shares use its separately snapshotted ISR interval of 998367 us.

```
isr_wake       1150.2/f 0.328/1.755/17.993 us 3.23%
  isr_pack       143.8/f 7.788/8.520/11.160 us 1.96%
  isr_dma_submit 143.8/f 0.613/0.920/7.925 us 0.21%
```

- Mean pack 8.520 us/call versus submit 0.920 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.23%, corresponding to 60.481 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.167×; Peak live show render: 0.894×; Transition / non-hold regime: 0.722×; Show end: 0.421×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `is_timeline_step` — 89.3% of the frame, 55.539 ms/f.
2. `is_draw_shape` — 89.0% of the frame, 55.384 ms/f.
3. `is_mesh_scan` — 74.1% of the frame, 46.091 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
