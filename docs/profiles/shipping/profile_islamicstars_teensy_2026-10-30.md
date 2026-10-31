# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with the matched full-cycle command below).
Raw capture: `C:\work\temp\ripple-unify-study-20261030\islamicstars-final-e55cab111.log`.
Matched baseline: `C:\work\temp\ripple-unify-study-20261030\islamicstars-before.log`. This replaces the 2026-10-26 report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz; COM3, serial 6362410; COM4 held as loader guard; flywheel + DMA ISRs live |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` mesh/transform hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars 288×144; source `e55cab1116daea41a4a9de5b793282040db0ed19` |
| Method | 210 s, window 16, `HS_PROFILE_TRANS_SPEED=4`, `HS_PROFILE_EPOCH_REVS=1920`; all 16 shapes and wrap |
| Reproduce | `HS_PROFILE_TREE="$PWD" bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` under the same pinned bench adapter |

Profile image size: `FLASH code/data/headers/free: 205256/197692/8700/1619968; RAM1 variables/code/padding/free: 313664/30808/1960/177856; RAM2 variables/free: 520064/4224`.

Full Phantasm image size: `FLASH code/data/headers/free: 1037344/749240/8488/236544; RAM1 variables/code/padding/free: 312736/196504/104/14944; RAM2 variables/free: 520064/4224`.

Exactness cross-check: frames 1473–1488, root 587564085 cycles / 600 = 979273.475 us versus wall sum 979274 us: **0.54 ppm**, strict canonical 5 ppm gate PASS.

All-window diagnostic: 0/208 windows exceed 5 ppm; maximum 4.74 ppm. This is reported separately from the canonical richest-window verdict.

Profile ELF SHA-256 `6cabb0c486a2095a1f4d7a3a6ed7dfa9d6c5bd7933062635a5cf29831293cf25`; Phantasm ELF SHA-256 `16f5d1c48cc7fecd69e1fd94a5a89a229d2f39f16201980dcb4d89f49749bad9`. The generator verifies retained ELF and compiler-envdump hashes against provenance.

## Frame cadence

**Baseline:** all 3328 complete frames: peak 54.046 ms at f249, mean 25.029 ms, spilled 0/3328. Ongoing frames 2–3328: peak 54.046 ms, mean 25.031 ms, spilled 0/3327. Startup f1: 18.674 ms. Excluded trailing raw rows: 9.

**Candidate:** all 3328 complete frames: peak 53.769 ms at f249, mean 25.023 ms, spilled 0/3328. Ongoing frames 2–3328: peak 53.769 ms, mean 25.025 ms, spilled 0/3327. Startup f1: 18.738 ms. Excluded trailing raw rows: 9.

Matched complete prefix: 3328 frames, ownership mismatches: 0. Ongoing peak delta: -0.277 ms; mean delta: -0.006 ms. The comparison includes construction, transitions, held geometry and wrap; it is not a hold-only result.

README cells: peak 🟢 53.769 (16), spilled 🟢 0/3327 (0.00%).
These cells use ongoing frames 2–3328; setup remains separate.

Observed ongoing peak margin: 8.731 ms against one display window. Required peak speedup: 1.000x (1.000x means no speedup required).

A display window is 62.5 ms at 480 RPM. Segment 0 renders an alternating 144×72 quadrant (10,368 pixels); the initial prime has a wider footprint. `canvas_buffer_wait` is display-sync idle. Render telemetry includes ISR time and excludes this wait.

## Phase-by-phase readout

Construction legs alternate with held solids; TS4 compresses the schedule. Representative windows maximize `is_timeline_step` within each stated regime. Trees preserve every logged row and nesting, including duplicate names. Percentages use the root frame counter; inclusive or tagged costs must not be summed.

### Startup (frames 1–16)

```text
frame                         57.133 ms 34.28 Mcyc 100.0%
  pov_preserve_half           137.47 us 82.48 kcyc   0.2%
    x0.9 146.64 us/call
  is_timeline_step            18.641 ms 11.18 Mcyc  32.6%
    hk_conway_previous         48.54 us 29.13 kcyc   0.1%
      x0.2 258.90 us/call
    is_build_draw             14.846 ms  8.91 Mcyc  26.0%
      is_build_scan           14.513 ms  8.71 Mcyc  25.4%
        scan_face_setup       112.70 us 67.62 kcyc   0.2% x1.8 64.40 us/call
        mesh_copy_capture      67.55 us 40.53 kcyc   0.1% x1.9 36.03 us/call
        mesh_capture          196.34 us 117.81 kcyc   0.3%
          x2.0 98.17 us/call
        scan_face_sweep        2.855 ms  1.71 Mcyc   5.0%
          scan_sweep_column   582.83 us 349.70 kcyc   1.0%
            x109.8 5.31 us/call
        mesh_sweep_paths       36.21 us 21.72 kcyc   0.1% x1.1 32.18 us/call
        mesh_classify         131.28 us 78.77 kcyc   0.2%
          x0.6 210.05 us/call
      is_mesh_transform        24.72 us 14.83 kcyc   0.0% x0.8 32.96 us/call
    hk_conway_compile          60.06 us 36.04 kcyc   0.1% x0.8 80.08 us/call
    hk_conway_sweep           118.52 us 71.11 kcyc   0.2%
      x0.8 158.03 us/call
    is_draw_shape              3.099 ms  1.86 Mcyc   5.4%
      is_mesh_scan             3.022 ms  1.81 Mcyc   5.3%
        scan_mesh_raster      12.697 ms  7.62 Mcyc  22.2%
          x23.5 540.29 us/call
        scan_face_setup       578.06 us 346.84 kcyc   1.0%
          x22.6 25.55 us/call
      is_face_offsets           4.30 us  2.58 kcyc   0.0% x0.2 17.22 us/call
      is_mesh_transform         6.51 us  3.90 kcyc   0.0% x0.2 26.03 us/call
      mesh_visibility          67.63 us 40.58 kcyc   0.1% x1.0 67.63 us/call
  is_ripple_prepare             3.70 us  2.22 kcyc   0.0% x1.0 3.70 us/call
  canvas_clear                100.57 us 60.34 kcyc   0.2%
    x1.0 100.57 us/call
  canvas_buffer_wait          38.204 ms 22.92 Mcyc  66.9%
    x1.0 38204.48 us/call
```

Wall min/avg/max=11.160/57.133/69.600 ms. Render mean/peak=18.929/23.782 ms. Display wait explains the gap between wall and render.

Tagged rows in this window: `scan_face_setup` (DUPLICATE-NAME/MIXED-PARENT), `is_mesh_transform` (DUPLICATE-NAME), `scan_mesh_raster` (MIXED-PARENT), `mesh_visibility` (MIXED-PARENT).

### Held solid (frames 513–528)

```text
scan_face_setup                96.20 us 57.72 kcyc   0.2% x1.9 51.31 us/call
mesh_copy_capture              62.30 us 37.38 kcyc   0.1% x2.6 23.73 us/call
mesh_capture                  567.57 us 340.54 kcyc   0.9%
  x1.6 363.25 us/call
scan_face_sweep                1.496 ms 897.86 kcyc   2.4%
  scan_sweep_column           538.57 us 323.14 kcyc   0.9%
    x40.3 13.36 us/call
mesh_sweep_paths               21.45 us 12.87 kcyc   0.0% x0.3 68.63 us/call
mesh_classify                 318.10 us 190.86 kcyc   0.5%
  x0.3 1017.93 us/call
frame                         62.328 ms 37.40 Mcyc 100.0%
  pov_preserve_half           145.59 us 87.35 kcyc   0.2%
    x1.0 145.59 us/call
  is_timeline_step            40.811 ms 24.49 Mcyc  65.5%
    is_draw_shape             40.655 ms 24.39 Mcyc  65.2%
      is_mesh_scan            36.187 ms 21.71 Mcyc  58.1%
        scan_mesh_raster      26.446 ms 15.87 Mcyc  42.4%
          x375.5 70.43 us/call
        scan_face_setup        5.033 ms  3.02 Mcyc   8.1%
          x385.1 13.07 us/call
      is_face_offsets          28.50 us 17.10 kcyc   0.0% x1.0 28.50 us/call
      is_mesh_transform        1.987 ms  1.19 Mcyc   3.2%
        x1.0 1986.58 us/call
      mesh_visibility          1.216 ms 729.38 kcyc   2.0%
        x1.0 1215.63 us/call
  is_ripple_prepare            17.01 us 10.21 kcyc   0.0% x1.0 17.01 us/call
  canvas_clear                 94.26 us 56.56 kcyc   0.2% x1.0 94.26 us/call
  canvas_buffer_wait          21.216 ms 12.73 Mcyc  34.0%
    x1.0 21215.97 us/call
```

Wall min/avg/max=47.421/62.328/72.582 ms. Render mean/peak=41.112/53.578 ms. Display wait explains the gap between wall and render.

Tagged rows in this window: `scan_face_setup` (DUPLICATE-NAME/MIXED-PARENT), `mesh_copy_capture` (MIXED-PARENT), `mesh_capture` (MIXED-PARENT), `scan_face_sweep` (MIXED-PARENT), `mesh_sweep_paths` (MIXED-PARENT), `mesh_classify` (MIXED-PARENT).

### Construction / transition (frames 1473–1488)

```text
scan_mesh_raster              27.072 ms 16.24 Mcyc  44.2%
  x271.7 99.64 us/call
scan_face_setup                2.837 ms  1.70 Mcyc   4.6%
  x278.4 10.19 us/call
mesh_visibility               558.91 us 335.35 kcyc   0.9%
  x1.0 558.91 us/call
frame                         61.205 ms 36.72 Mcyc 100.0%
  pov_preserve_half           141.93 us 85.16 kcyc   0.2%
    x1.0 141.93 us/call
  is_timeline_step            41.626 ms 24.98 Mcyc  68.0%
    hk_conway_previous        629.88 us 377.93 kcyc   1.0%
      x0.4 1439.72 us/call
    is_build_draw             38.783 ms 23.27 Mcyc  63.4%
      is_build_scan           37.384 ms 22.43 Mcyc  61.1%
        scan_face_setup       655.37 us 393.22 kcyc   1.1%
          x15.9 41.28 us/call
        mesh_copy_capture     303.97 us 182.38 kcyc   0.5%
          x18.2 16.66 us/call
        mesh_capture          747.79 us 448.67 kcyc   1.2%
          x3.6 206.29 us/call
        scan_face_sweep        1.693 ms  1.02 Mcyc   2.8%
          scan_sweep_column   529.55 us 317.73 kcyc   0.9%
            x121.9 4.34 us/call
        mesh_sweep_paths       31.76 us 19.06 kcyc   0.1% x1.4 22.10 us/call
        mesh_classify         496.41 us 297.85 kcyc   0.8%
          x0.8 661.88 us/call
      is_mesh_transform        74.30 us 44.58 kcyc   0.1% x1.0 74.30 us/call
    hk_conway_compile         379.25 us 227.55 kcyc   0.6%
      x1.0 379.25 us/call
    hk_conway_sweep           891.77 us 535.06 kcyc   1.5%
      x1.0 891.77 us/call
  is_ripple_prepare             3.62 us  2.17 kcyc   0.0% x1.0 3.62 us/call
  canvas_clear                 94.77 us 56.86 kcyc   0.2% x1.0 94.77 us/call
  canvas_buffer_wait          19.111 ms 11.47 Mcyc  31.2%
    x1.0 19111.50 us/call
```

Wall min/avg/max=39.632/61.204/73.425 ms. Render mean/peak=42.093/50.410 ms. Display wait explains the gap between wall and render.

Tagged rows in this window: `scan_mesh_raster` (MIXED-PARENT), `scan_face_setup` (DUPLICATE-NAME/MIXED-PARENT), `mesh_visibility` (MIXED-PARENT).

### Per-preset table

All 16 shapes are observed; wrap to the first occurs at f1224. Rows are ranked by the chosen held-window `is_timeline_step` cost. **Spawn V/E/F/I are seed geometry, not the constructed mesh.** Final `Built Shape` geometry is retained separately in JSON.

The primary cadence columns use every complete frame owned by each shape, including construction and transitions. Ownership begins with the frame after its `Spawning Shape` marker. A boundary window may contain two owners; its individual frames remain assigned correctly.

| # | Shape | Seed V/E/F/I | Frames | Peak render ms | Spilled | Construction / pure-owner windows |
|---|---|---|--:|--:|--:|--:|
| 6 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 20/30/12/60 | 264 | 53.578 | 0/264 | 7/13 |
| 4 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 60/90/32/180 | 216 | 44.891 | 0/216 | 6/10 |
| 15 | icosahedron_truncate25_hk60_ambo_relax100_hk60 | 12/30/20/60 | 160 | 41.562 | 0/160 | 5/8 |
| 8 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 120/180/62/360 | 234 | 51.371 | 0/234 | 6/10 |
| 3 | truncatedCube_hk90_meta | 24/36/14/72 | 276 | 53.769 | 0/276 | 9/14 |
| 14 | icosahedron_snub_relax_truncate033_hankin62 | 12/30/20/60 | 144 | 38.053 | 0/144 | 4/7 |
| 1 | truncatedIcosahedron_hk58_chamfer63 | 60/90/32/180 | 186 | 34.836 | 0/186 | 3/9 |
| 13 | dodecahedron_hk72_ambo_dual_hk20 | 20/30/12/60 | 176 | 35.022 | 0/176 | 6/9 |
| 11 | rhombicuboctahedron_hk63_ambo_hk63 | 24/48/26/96 | 177 | 36.829 | 0/177 | 5/9 |
| 12 | dodecahedron_hk54_ambo_hk72 | 20/30/12/60 | 140 | 31.034 | 0/140 | 3/6 |
| 0 | dodecahedron_hk62_ambo_hk62 | 20/30/12/60 | 209 | 31.306 | 0/209 | 6/11 |
| 7 | octahedron_hk17_ambo_hk73 | 6/12/8/24 | 210 | 36.688 | 0/210 | 5/11 |
| 10 | snubDodecahedron_truncate5d_ambo_dual | 60/150/92/300 | 234 | 36.242 | 0/234 | 6/12 |
| 5 | icosahedron_ambo_truncate033_hankin59 | 12/30/20/60 | 204 | 35.809 | 0/204 | 5/11 |
| 2 | dodecahedron_bevel2_relax_gyro | 20/30/12/60 | 264 | 34.873 | 0/264 | 7/13 |
| 9 | icosidodecahedron_truncate5d_ambo_dual | 30/60/32/120 | 234 | 35.055 | 0/234 | 8/12 |

Matched-prefix per-shape comparison (all owned frames, not only held geometry):

| Shape | Before / after frames | Before peak ms | After peak ms | Delta ms | Before / after spilled |
|---|--:|--:|--:|--:|--:|
| dodecahedron_hk62_ambo_hk62 | 209/209 | 31.146 | 31.306 | +0.160 | 0/0 |
| truncatedIcosahedron_hk58_chamfer63 | 186/186 | 34.806 | 34.836 | +0.030 | 0/0 |
| dodecahedron_bevel2_relax_gyro | 264/264 | 34.869 | 34.873 | +0.004 | 0/0 |
| truncatedCube_hk90_meta | 276/276 | 54.046 | 53.769 | -0.277 | 0/0 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 216/216 | 44.529 | 44.891 | +0.362 | 0/0 |
| icosahedron_ambo_truncate033_hankin59 | 204/204 | 35.795 | 35.809 | +0.014 | 0/0 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 264/264 | 53.602 | 53.578 | -0.024 | 0/0 |
| octahedron_hk17_ambo_hk73 | 210/210 | 36.741 | 36.688 | -0.053 | 0/0 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 234/234 | 51.807 | 51.371 | -0.436 | 0/0 |
| icosidodecahedron_truncate5d_ambo_dual | 234/234 | 35.041 | 35.055 | +0.014 | 0/0 |
| snubDodecahedron_truncate5d_ambo_dual | 234/234 | 36.394 | 36.242 | -0.152 | 0/0 |
| rhombicuboctahedron_hk63_ambo_hk63 | 177/177 | 36.685 | 36.829 | +0.144 | 0/0 |
| dodecahedron_hk54_ambo_hk72 | 140/140 | 30.863 | 31.034 | +0.171 | 0/0 |
| dodecahedron_hk72_ambo_dual_hk20 | 176/176 | 34.763 | 35.022 | +0.259 | 0/0 |
| icosahedron_snub_relax_truncate033_hankin62 | 144/144 | 37.927 | 38.053 | +0.126 | 0/0 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 160/160 | 41.323 | 41.562 | +0.239 | 0/0 |

Held-window detail: one owner, no `is_build_draw`, `is_draw_shape` calls equal the 16 window frames, then modal `scan_mesh_raster` calls within that shape. No fallback to construction windows is allowed. The highest-cost remaining window supplies this table.

| Shape | Window | Clean / held windows | Blended px/f | Dominant ms | Render ms | Observed fps |
|---|---|--:|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 497–512 | 1/6 | unavailable | 34.790 | 35.082 | 16.06 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 353–368 | 1/4 | unavailable | 34.406 | 34.702 | 15.85 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 1185–1200 | 1/3 | unavailable | 28.857 | 29.154 | 15.95 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 657–672 | 1/4 | unavailable | 27.915 | 28.210 | 15.96 |
| truncatedCube_hk90_meta | 273–288 | 1/5 | unavailable | 27.009 | 27.310 | 16.03 |
| icosahedron_snub_relax_truncate033_hankin62 | 1105–1120 | 1/3 | unavailable | 25.487 | 25.786 | 15.98 |
| truncatedIcosahedron_hk58_chamfer63 | 97–112 | 1/6 | unavailable | 24.950 | 25.243 | 16.03 |
| dodecahedron_hk72_ambo_dual_hk20 | 1041–1056 | 1/3 | unavailable | 23.945 | 24.240 | 15.92 |
| rhombicuboctahedron_hk63_ambo_hk63 | 881–896 | 1/4 | unavailable | 23.733 | 24.028 | 16.01 |
| dodecahedron_hk54_ambo_hk72 | 945–960 | 1/3 | unavailable | 22.327 | 22.625 | 15.99 |
| dodecahedron_hk62_ambo_hk62 | 33–48 | 1/5 | unavailable | 21.563 | 21.860 | 16.02 |
| octahedron_hk17_ambo_hk73 | 1809–1824 | 2/6 | unavailable | 20.521 | 20.827 | 16.05 |
| snubDodecahedron_truncate5d_ambo_dual | 801–816 | 1/6 | unavailable | 20.067 | 20.363 | 16.21 |
| icosahedron_ambo_truncate033_hankin59 | 417–432 | 1/6 | unavailable | 19.819 | 20.119 | 16.00 |
| dodecahedron_bevel2_relax_gyro | 177–192 | 1/6 | unavailable | 19.337 | 19.631 | 16.23 |
| icosidodecahedron_truncate5d_ambo_dual | 737–752 | 1/4 | unavailable | 16.770 | 17.071 | 15.98 |

### Per-pixel figures

Timing images do not contain `filter_blend` counters. Blended-pixel coverage, cycles/blend and scan cycles/blend are unavailable; the 10,368-pixel quadrant is a geometric area, not a measured blend count.

## Column-ISR / DMA marshaling cost

Aggregated across all 208 complete ISR windows (208.401 s). ISR intervals include their own boundary/readout elapsed time.

```text
isr_wake        1154.30/f  18433.2/s
  min/avg/max 0.455/1.980/32.517 us; CPU 3.65%
isr_pack         143.93/f   2298.5/s
  min/avg/max 5.992/9.318/14.732 us; CPU 2.14%
isr_dma_submit   143.93/f   2298.5/s
  min/avg/max 0.578/0.932/9.043 us; CPU 0.21%
```

- Submit is CPU setup time; pack performs the larger measured CPU work. Async wire-transfer duration is not measured.
- Inclusive `isr_wake` uses 3.65% of the ISR interval, leaving an estimated 60.218 ms outside measured wake work per display window. Pack and submit are nested; their shares are not added to wake.
- Render telemetry already includes interrupts; cadence uses the 62.5 ms window directly. DMA-completion and other interrupts are not separately measured.

## Summary ranking

1. `is_timeline_step` — 39.6% of root cycles, 24.684 ms/frame inclusive.
2. `is_draw_shape` — 23.1% of root cycles, 14.424 ms/frame inclusive.
3. `is_mesh_scan` — 22.0% of root cycles, 13.747 ms/frame inclusive.

These inclusive scopes overlap. Native captures establish unchanged pixel output; their timings are not substituted for this device capture.

## Caveats

- CYCCNT scopes include ISR time.
- `MIXED-PARENT` and `DUPLICATE-NAME` rows preserve logged placement, but do not identify exclusive work. No counter rows are discarded from the JSON trees.
- `filter_blend` parenting can place its subtree under the first caller; these timing captures have no blend counters.
- Scopes are per stage/call rather than newly added per-pixel probes. Instrumentation has overhead.
- The image uses shipping selective-O3 regions in mesh scan and spherical transforms, with an `-Os` base.
- TS4 compresses dwell/construction scheduling; it is fixed in both captures and is not representative of the default-duration cycle.
- Source SHA, retained ELF hashes and envdump hashes are recorded per capture. The generator does not infer capture-time state from a later checkout.
- Complete counter windows define the report span; trailing raw telemetry is recorded in JSON and excluded consistently.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`, `HS_PROFILE_TRANS_SPEED=4`, `HS_PROFILE_EPOCH_REVS=1920`. `just profile IslamicStars` invokes the supported capture flow; use the full-cycle 210 s/window 16 command above for this comparison.
