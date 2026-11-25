# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-11-17, **global -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-islamic-o3-evidence-20261117\islamicstars_o3.log`. Supersedes the [historical TS4 capture](profile_islamicstars_historical_ts4_teensy_2026-09-28.md) as the ranked O3 result.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile_o3`; global `-O3`, newlib-nano; same `-ffast-math -fno-finite-math-only` base flags |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 306 s / 4896 nominal frames at 16 fps; window 16. Raw capture 309 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh IslamicStars profile_o3 309 16`; analyze only the 306 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:213928, data:196476, headers:8412; free:1612800 / RAM1: variables:314720, code:40136, padding:25400; free:144032 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 2049–2064 root cycles / 600 MHz = 1002747.088 us; measured wall sum 1002749 us, delta 1.91 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 10:34, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `3ca27a4713a68867160a20b33750076f13b72f20439bd8b2238591fe83d1457e`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 306 s (4896 nominal frames at 16 fps). Show-frame peak render is **57.095 ms**, spilled **0/4896**. Setup frame 1 is **16.469 ms** and remains separate. Live show peak is **57.095 ms**, spilled **0/4895**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–4896. The raw log retains 4921 frame records; the 25 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 57.095 (16), spilled 🟢 0/4895 (0.00%).

README setup frame 1: 16.469 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    56.308 ms  33.785 Mcyc 100.0%
  pov_preserve_half      137.00 us   82.20 kcyc   0.2% x0.9375 146 us/call
  is_timeline_step       10.231 ms   6.139 Mcyc  18.2%
    is_draw_shape        10.037 ms   6.022 Mcyc  17.8%
      is_mesh_scan        9.730 ms   5.838 Mcyc  17.3%
        scan_mesh_raster   9.373 ms   5.624 Mcyc  16.6% x7.9375 1181 us/call
        scan_face_setup  157.53 us   94.52 kcyc   0.3% x7.9375 20 us/call
      is_face_offsets     42.35 us   25.41 kcyc   0.1% x1 42 us/call
      is_mesh_transform   37.75 us   22.65 kcyc   0.1% x1 38 us/call
      mesh_visibility     45.24 us   27.14 kcyc   0.1% x1 45 us/call
  is_ripple_prepare        1.53 us    0.92 kcyc   0.0% x1 2 us/call
  canvas_clear            98.66 us   59.20 kcyc   0.2% x1 99 us/call
  canvas_buffer_wait     45.801 ms  27.481 Mcyc  81.3% x1 45801 us/call
```

Wall min/avg/max = 9.426/56.307/63.772 ms. `is_timeline_step` averages 10.231 ms/f; total render averages 10.506 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 2049–2064)

```
scan_face_setup          239.79 us  143.87 kcyc   0.4% x12.5 19 us/call
mesh_copy_capture        180.84 us  108.51 kcyc   0.3% x16.375 11 us/call
mesh_capture              2.879 ms   1.727 Mcyc   4.6% x4.875 591 us/call
mesh_classify             1.447 ms  868.32 kcyc   2.3% x1 1447 us/call
frame                    62.672 ms  37.603 Mcyc 100.0%
  pov_preserve_half      140.72 us   84.43 kcyc   0.2% x1 141 us/call
  is_timeline_step       54.294 ms  32.576 Mcyc  86.6%
    is_draw_shape        54.122 ms  32.473 Mcyc  86.4%
      is_mesh_scan       44.951 ms  26.971 Mcyc  71.7%
        scan_face_sweep   3.827 ms   2.296 Mcyc   6.1%
          scan_sweep_column   1.382 ms  829.33 kcyc   2.2% x90 15 us/call
        mesh_sweep_paths  123.34 us   74.01 kcyc   0.2% x0.6875 179 us/call
        scan_mesh_raster  27.599 ms  16.560 Mcyc  44.0% x393.812 70 us/call
        scan_face_setup   5.569 ms   3.342 Mcyc   8.9% x436.562 13 us/call
      is_face_offsets     59.11 us   35.47 kcyc   0.1% x1 59 us/call
      is_mesh_transform   2.418 ms   1.451 Mcyc   3.9% x1 2418 us/call
      mesh_visibility     1.225 ms  735.10 kcyc   2.0% x1 1225 us/call
  is_ripple_prepare       18.32 us   10.99 kcyc   0.0% x1 18 us/call
  canvas_clear            95.83 us   57.50 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait      8.054 ms   4.832 Mcyc  12.9% x1 8054 us/call
```

Wall min/avg/max = 60.052/62.671/65.646 ms. `is_timeline_step` averages 54.294 ms/f; total render averages 54.618 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 2481–2496)

```
scan_mesh_raster         33.493 ms  20.096 Mcyc  53.4% x197.812 169 us/call
scan_face_setup           2.702 ms   1.621 Mcyc   4.3% x207.75 13 us/call
mesh_visibility          533.80 us  320.28 kcyc   0.9% x1 534 us/call
frame                    62.774 ms  37.664 Mcyc 100.0%
  pov_preserve_half      140.66 us   84.39 kcyc   0.2% x1 141 us/call
  is_timeline_step       42.736 ms  25.642 Mcyc  68.1%
    hk_conway_previous   860.00 us  516.00 kcyc   1.4% x0.9375 917 us/call
    is_build_draw        40.266 ms  24.160 Mcyc  64.1%
      is_build_scan      39.039 ms  23.423 Mcyc  62.2%
        scan_face_setup  128.85 us   77.31 kcyc   0.2% x7.125 18 us/call
        mesh_copy_capture   26.15 us   15.69 kcyc   0.0% x7.125 4 us/call
        mesh_capture     277.29 us  166.37 kcyc   0.4% x2.8125 99 us/call
        mesh_classify    326.33 us  195.80 kcyc   0.5% x0.5625 580 us/call
      is_mesh_transform   74.42 us   44.65 kcyc   0.1% x1 74 us/call
    hk_conway_compile    364.62 us  218.77 kcyc   0.6% x1 365 us/call
    hk_conway_sweep      985.61 us  591.36 kcyc   1.6% x1 986 us/call
  is_ripple_prepare        3.05 us    1.83 kcyc   0.0% x1 3 us/call
  canvas_clear            92.17 us   55.30 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     19.211 ms  11.527 Mcyc  30.6% x1 19211 us/call
```

Wall min/avg/max = 54.675/62.773/73.138 ms. `is_timeline_step` averages 42.736 ms/f; total render averages 43.563 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 4881–4896)

```
frame                    62.410 ms  37.446 Mcyc 100.0%
  pov_preserve_half      146.19 us   87.71 kcyc   0.2% x1 146 us/call
  is_timeline_step       25.946 ms  15.567 Mcyc  41.6%
    is_draw_shape        25.713 ms  15.428 Mcyc  41.2%
      is_mesh_scan       24.230 ms  14.538 Mcyc  38.8%
        scan_mesh_raster  21.021 ms  12.613 Mcyc  33.7% x178.125 118 us/call
        scan_face_setup   2.412 ms   1.447 Mcyc   3.9% x187.938 13 us/call
      is_face_offsets    307.18 us  184.31 kcyc   0.5% x1 307 us/call
      is_mesh_transform  106.83 us   64.10 kcyc   0.2% x1 107 us/call
      mesh_visibility    599.99 us  359.99 kcyc   1.0% x1 600 us/call
  is_ripple_prepare        1.60 us    0.96 kcyc   0.0% x1 2 us/call
  canvas_clear            88.17 us   52.90 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     36.185 ms  21.711 Mcyc  58.0% x1 36185 us/call
```

Wall min/avg/max = 60.730/62.409/63.853 ms. `is_timeline_step` averages 25.946 ms/f; total render averages 26.224 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Qualified windows have one shape owner, one is_draw_shape call per frame, no is_build_draw calls, and lie after the source-derived build/fade-in interval and before fade-out. The current visibility pass culls faces before scan_mesh_raster, so that counter varies with orientation and is not equal to the complete face count. Seed and completed geometry are listed separately. Completed geometry remains animated by the shipping ripples and orientation.

The show-scoped data covers 16 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 9 | 2049–2064 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 3240/4320/1082/8640 | — | 54.294 | 54.618 | 16 | 57.095 | 0/352 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 9 | 1425–1440 | built-shape full-envelope live-ripple window | seed 60/90/32/180; built 1620/2160/542/4320 | — | 46.256 | 46.573 | 16 | 48.107 | 0/288 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 9 | 4785–4800 | built-shape full-envelope live-ripple window | seed 12/30/20/60; built 1620/2160/542/4320 | — | 38.231 | 38.546 | 16 | 41.084 | 0/320 |
| icosahedron_snub_relax_truncate033_hankin62 | 9 | 4465–4480 | built-shape full-envelope live-ripple window | seed 12/30/20/60; built 1350/1800/452/3600 | — | 34.587 | 34.902 | 16 | 37.526 | 0/288 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 9 | 2641–2656 | built-shape full-envelope live-ripple window | seed 120/180/62/360; built 542/1080/540/2160 | — | 33.626 | 33.945 | 16 | 47.681 | 0/312 |
| truncatedCube_hk90_meta | 9 | 1169–1184 | built-shape full-envelope live-ripple window | seed 24/36/14/72; built 218/648/432/1296 | — | 32.490 | 32.797 | 16 | 48.648 | 0/368 |
| truncatedIcosahedron_hk58_chamfer63 | 9 | 417–432 | built-shape full-envelope live-ripple window | seed 60/90/32/180; built 990/1440/452/2880 | — | 30.928 | 31.243 | 16 | 32.855 | 0/248 |
| dodecahedron_hk72_ambo_dual_hk20 | 9 | 4177–4192 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 540/720/182/1440 | — | 30.759 | 31.072 | 16 | 33.362 | 0/352 |
| rhombicuboctahedron_hk63_ambo_hk63 | 9 | 3537–3552 | built-shape full-envelope live-ripple window | seed 24/48/26/96; built 864/1152/290/2304 | — | 29.639 | 29.958 | 16 | 32.651 | 0/280 |
| dodecahedron_hk54_ambo_hk72 | 9 | 3825–3840 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 540/720/182/1440 | — | 29.588 | 29.907 | 16 | 32.454 | 0/280 |
| icosahedron_ambo_truncate033_hankin59 | 9 | 1697–1712 | built-shape full-envelope live-ripple window | seed 12/30/20/60; built 540/720/182/1440 | — | 27.720 | 28.052 | 16 | 32.782 | 0/272 |
| snubDodecahedron_truncate5d_ambo_dual | 9 | 3265–3280 | built-shape full-envelope live-ripple window | seed 60/150/92/300; built 452/900/450/1800 | — | 25.925 | 26.242 | 16 | 33.751 | 0/312 |
| dodecahedron_hk62_ambo_hk62 | 9 | 177–192 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 540/720/182/1440 | — | 25.029 | 25.343 | 16 | 26.069 | 0/279 |
| dodecahedron_bevel2_relax_gyro | 9 | 769–784 | built-shape full-envelope live-ripple window | seed 20/30/12/60; built 542/900/360/1800 | — | 24.592 | 24.913 | 16 | 33.376 | 0/352 |
| octahedron_hk17_ambo_hk73 | 9 | 2321–2336 | built-shape full-envelope live-ripple window | seed 6/12/8/24; built 216/288/74/576 | — | 24.362 | 24.680 | 16 | 33.471 | 0/280 |
| icosidodecahedron_truncate5d_ambo_dual | 9 | 2945–2960 | built-shape full-envelope live-ripple window | seed 30/60/32/120; built 182/360/180/720 | — | 18.788 | 19.103 | 16 | 26.766 | 0/312 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 2049–2064; CPU shares use its separately snapshotted ISR interval of 1005664 us.

```
isr_wake       1158.6/f 0.225/1.621/22.385 us 2.99%
  isr_pack       144.8/f 7.831/8.498/20.843 us 1.95%
  isr_dma_submit 144.8/f 0.596/0.918/7.595 us 0.21%
```

- Mean pack 8.498 us/call versus submit 0.918 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.99%, corresponding to 60.631 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.168×; Peak live show render: 0.874×; Transition / non-hold regime: 0.697×; Show end: 0.420×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `is_timeline_step` — 86.6% of the frame, 54.294 ms/f.
2. `is_draw_shape` — 86.4% of the frame, 54.122 ms/f.
3. `is_mesh_scan` — 71.7% of the frame, 44.951 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Global O3 covers the translation units; the shipping-only `HS_O3` annotations compile out. This single-effect reference is not the full-roster shipping image.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.

## Matched shipping comparison

Same clean source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21`, Teensy COM3 (serial 6362410), default Trans Speed, 16 recipes, 16-frame windows, and the same 306-second authored show budget (frames 2–4896 for live statistics). The shipping capture was taken at 09:09 on the same day. Each configuration has one capture.

| Metric | Shipping selective O3 | Global O3 |
|---|--:|--:|
| Live peak ms | 58.059 | 57.095 |
| Live mean ms | 24.548 | 24.699 |
| Live spilled | 0/4895 | 0/4895 |

Peak improves by 0.964 ms (1.66%). Both peaks occur on the 1,082-face dodecahedron recipe. The historical September 28 result used TS4, 23 recipes and older rendering code; its 49.112 ms peak is not a matched optimization delta against current shipping.

The following costs average the same nine qualified show windows of that heavy recipe. Rows are inclusive and must not be summed.

| Scope | Shipping ms/f | O3 ms/f | Saved ms/f |
|---|--:|--:|--:|
| `is_timeline_step` | 43.432 | 42.700 | +0.731 |
| `is_mesh_scan` | 37.906 | 37.322 | +0.584 |
| `scan_mesh_raster` | 26.706 | 26.667 | +0.039 |
| `scan_face_setup` | 5.188 | 4.986 | +0.203 |
| `scan_face_sweep` | 1.733 | 1.531 | +0.202 |
| `scan_sweep_column` | 0.706 | 0.571 | +0.135 |
| `mesh_capture` | 1.126 | 1.185 | -0.059 |

In the shipping peak window (2049–2064), swept-face scanning falls from 5.529 to 3.827 ms/f, but its calls also fall from 1.000 to 0.688 per frame; sweep-column calls fall from 143.25 to 90.00 per frame. This is a change in path mix as well as timing, not proof of a 1.7 ms core compiler win. The main face raster costs 26.916 versus 27.599 ms/f in that window, and 26.706 versus 26.667 ms/f across the nine qualified heavy windows.

Other per-recipe peaks move in both directions: truncatedCube_hk90_meta improves from 52.982 to 48.648 ms, while dodecahedron_hk54_ambo_hk72 rises from 26.867 to 32.454 ms. These single-capture maxima do not establish a uniformly faster core path.

The O3 profile image costs +28,528 B FLASH code and +9,456 B ITCM code; RAM1 free space changes by -32,800 B because the image crosses a 32 KiB ITCM allocation boundary. These are single-effect image deltas, not a measured full-roster selective change.

Core raster functions already receive selective O3. The shipping attributes disable loop unswitching to control code growth; global O3 restores it and changes inlining/layout. The capture does not demonstrate a large core saving worth importing. Any narrow unswitching experiment would need its own peak, render-equivalence and full-image memory qualification.

[Current shipping report](../shipping/profile_islamicstars_teensy_2026-11-20.md). [Historical TS4 capture](profile_islamicstars_historical_ts4_teensy_2026-09-28.md).
