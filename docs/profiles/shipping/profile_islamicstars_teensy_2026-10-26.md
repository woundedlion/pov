# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-26, **selective -O3**)

Point-in-time snapshot (regenerate with the full-cycle TS4 command in [Setup](#setup)).
Raw capture: `build/prof/islamicstars_review45_after.log`, captured 2026-10-26 09:28 America/Los_Angeles on COM4.
Matched baseline: `build/prof/islamicstars_review45_before.log`. Replaces `profile_islamicstars_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel and DMA ISRs live |
| Image | `profile`, shipping `-Os` base; selective `HS_O3` mesh scan regions and IslamicStars transforms |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, source `771290aa484e3b99b30b049030ab145d61ade1d9` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 210 s; `HS_PROFILE_TRANS_SPEED=4`, `HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_EPOCH_REVS=1920 -D HS_PROFILE_TRANS_SPEED=4"` |

Profile image size: `FLASH: code:190424, data:208884, headers:8244, free:1624064 / RAM1: variables:313664, code:30840, padding:1928, free:177856 / RAM2: variables:520064, free:4224`.

Full Phantasm image size: `FLASH: code:983472, data:853320, headers:8456, free:186368 / RAM1: variables:312736, code:176040, padding:20568, free:14944 / RAM2: variables:520064, free:4224`.

Baseline source `d459cc1969e465d5c208b64165cdf8fa4425d2c4`. Profile ELF SHA-256 `ebd99596297c4c71af64cb36a7a5001d091b584a74a11e24ad0ac57a44dd3be6`; Phantasm ELF SHA-256 `a88266b6028854037bf1e3435cc8e18e48e3e7f636fb152e8aff54474a87628c`.

Adjacent provenance, profile/Phantasm build logs and compiler envdumps are copied beside both raw captures. Provenance retains the original ELF artifact paths.

The integrated code tip is `eac6c22a58cfaf313725c9c41960f50fde2101ed`. Its differences from the captured source are confined to `tests/test_death.h` and `tests/test_plot_parity.h`; firmware source is identical.

Exactness cross-check: frames 1473–1488, root 587584496 cycles ÷ 600 = 979307.493 µs versus measured wall sum 979542 µs: **239.4 ppm**. The strict 5 ppm gate fails; this report does not claim fully validated cycle-counter exactness.

The harness timestamps outside the inner `frame` RAII scope. Boundary overhead and interrupts between these boundaries are not isolated by the capture. No parser threshold was changed.

## Frame cadence

**Steady render, frame 1 excluded:** peak **53.854 ms**, mean **25.159 ms**, **0/3327 spilled frames** (0.00%).
Complete counter windows cover frames 1–3328; transitions remain included. Setup frame 1 renders 18.791 ms. 8 trailing raw frames are excluded from aggregate and ownership tables.

Pass aggregate including setup (`parse_profile.py ... windows`): peak 53.854 ms; 0/3328 spilled. `is_timeline_step` averages 24.785 ms/frame; its richest window averages 41.779 ms/frame (frames 1473–1488).

A display window is 62.5 ms at 480 RPM. Segment 0 renders one alternating 144×72 quadrant after the initial full-width half-height prime; nominal cadence is 16 fps when render stays within the window. `canvas_buffer_wait` is the round-up idle to the next display flip.

The observed steady peak leaves 8.646 ms against the display budget. No cadence speedup is required for the observed held-solid and construction phases.

## Phase-by-phase readout

Phase schedule: construction legs alternate with held solids; TS4 compresses their dwell. Inclusive scopes overlap; MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs. Trees show recorded nesting and per-frame averages; leaf call counts are reported separately.

### Startup (window frames 1–16)

```text
frame                              57.131 ms  34.28 Mcyc 100.0%
  x1.00  57130.512 us/call
  pov_preserve_half                 0.135 ms   0.08 Mcyc   0.2%
    x0.94  144.203 us/call
  is_timeline_step                 18.691 ms  11.21 Mcyc  32.7%
    x1.00  18691.481 us/call
    hk_conway_previous              0.048 ms   0.03 Mcyc   0.1%
      x0.19  258.412 us/call
    is_build_draw                  14.892 ms   8.94 Mcyc  26.1%
      x0.75  19855.800 us/call
      is_build_scan                14.544 ms   8.73 Mcyc  25.5%
        x0.75  19391.544 us/call
        scan_face_setup             0.082 ms   0.05 Mcyc   0.1%
          x1.75  46.900 us/call
        mesh_copy_capture           0.060 ms   0.04 Mcyc   0.1%
          x1.88  31.962 us/call
        mesh_capture                0.215 ms   0.13 Mcyc   0.4%
          x2.00  107.446 us/call
        scan_face_sweep             2.882 ms   1.73 Mcyc   5.0%
          x1.12  2561.455 us/call
          scan_sweep_column         0.594 ms   0.36 Mcyc   1.0%
            x109.75  5.416 us/call
        mesh_sweep_paths            0.044 ms   0.03 Mcyc   0.1%
          x1.12  39.121 us/call
        mesh_classify               0.124 ms   0.07 Mcyc   0.2%
          x0.62  198.292 us/call
      is_mesh_transform             0.026 ms   0.02 Mcyc   0.0%
        x0.75  34.339 us/call
    hk_conway_compile               0.052 ms   0.03 Mcyc   0.1%
      x0.75  68.957 us/call
    hk_conway_sweep                 0.119 ms   0.07 Mcyc   0.2%
      x0.75  158.463 us/call
    is_draw_shape                   3.117 ms   1.87 Mcyc   5.5%
      x0.25  12466.454 us/call
      is_mesh_scan                  3.041 ms   1.82 Mcyc   5.3%
        x0.25  12164.588 us/call
        scan_mesh_raster           12.724 ms   7.63 Mcyc  22.3%
          x23.50  541.439 us/call
        scan_face_setup             0.564 ms   0.34 Mcyc   1.0%
          x22.62  24.939 us/call
      is_face_offsets               0.002 ms   0.00 Mcyc   0.0%
        x0.25  9.209 us/call
      is_mesh_transform             0.006 ms   0.00 Mcyc   0.0%
        x0.25  22.907 us/call
      mesh_visibility               0.066 ms   0.04 Mcyc   0.1%
        x1.00  66.231 us/call
  is_ripple_prepare                 0.003 ms   0.00 Mcyc   0.0%
    x1.00  2.641 us/call
  canvas_clear                      0.101 ms   0.06 Mcyc   0.2%
    x1.00  101.381 us/call
  canvas_buffer_wait               38.135 ms  22.88 Mcyc  66.8%
    x1.00  38134.933 us/call
```

Wall min/avg/max = 11.219/57.144/69.440 ms. Render window mean is 18.996 ms; its peak frame is 23.793 ms. The difference from wall includes display-sync idle.

`is_timeline_step`: 1.00 calls/frame, 18691.481 µs/call. Scope tags, where present, retain the parser limitations.

### Held solid (window frames 513–528)

```text
scan_face_setup                     0.067 ms   0.04 Mcyc   0.1%
  x1.88  35.593 us/call
mesh_copy_capture                   0.052 ms   0.03 Mcyc   0.1%
  x2.62  19.703 us/call
mesh_capture                        0.565 ms   0.34 Mcyc   0.9%
  x1.56  361.439 us/call
scan_face_sweep                     1.486 ms   0.89 Mcyc   2.4%
  x0.31  4754.509 us/call
  scan_sweep_column                 0.536 ms   0.32 Mcyc   0.9%
    x40.31  13.304 us/call
mesh_sweep_paths                    0.023 ms   0.01 Mcyc   0.0%
  x0.31  73.824 us/call
mesh_classify                       0.324 ms   0.19 Mcyc   0.5%
  x0.31  1038.242 us/call
frame                              62.300 ms  37.38 Mcyc 100.0%
  x1.00  62300.187 us/call
  pov_preserve_half                 0.143 ms   0.09 Mcyc   0.2%
    x1.00  143.381 us/call
  is_timeline_step                 40.978 ms  24.59 Mcyc  65.8%
    x1.00  40978.433 us/call
    is_draw_shape                  40.801 ms  24.48 Mcyc  65.5%
      x1.00  40800.939 us/call
      is_mesh_scan                 36.268 ms  21.76 Mcyc  58.2%
        x1.00  36268.405 us/call
        scan_mesh_raster           26.560 ms  15.94 Mcyc  42.6%
          x375.50  70.733 us/call
        scan_face_setup             5.027 ms   3.02 Mcyc   8.1%
          x385.12  13.053 us/call
      is_face_offsets               0.020 ms   0.01 Mcyc   0.0%
        x1.00  20.104 us/call
      is_mesh_transform             2.062 ms   1.24 Mcyc   3.3%
        x1.00  2062.060 us/call
      mesh_visibility               1.195 ms   0.72 Mcyc   1.9%
        x1.00  1195.034 us/call
  is_ripple_prepare                 0.015 ms   0.01 Mcyc   0.0%
    x1.00  14.535 us/call
  canvas_clear                      0.095 ms   0.06 Mcyc   0.2%
    x1.00  95.091 us/call
  canvas_buffer_wait               21.010 ms  12.61 Mcyc  33.7%
    x1.00  21010.163 us/call
```

Wall min/avg/max = 47.436/62.313/72.394 ms. Render window mean is 41.290 ms; its peak frame is 53.665 ms. The difference from wall includes display-sync idle.

`is_timeline_step`: 1.00 calls/frame, 40978.433 µs/call. Scope tags, where present, retain the parser limitations.

### Construction / transition (window frames 1473–1488)

```text
scan_mesh_raster                   27.127 ms  16.28 Mcyc  44.3%
  x271.69  99.845 us/call
scan_face_setup                     2.774 ms   1.66 Mcyc   4.5%
  x278.44  9.963 us/call
mesh_visibility                     0.562 ms   0.34 Mcyc   0.9%
  x1.00  561.699 us/call
frame                              61.207 ms  36.72 Mcyc 100.0%
  x1.00  61206.718 us/call
  pov_preserve_half                 0.139 ms   0.08 Mcyc   0.2%
    x1.00  139.063 us/call
  is_timeline_step                 41.779 ms  25.07 Mcyc  68.3%
    x1.00  41778.812 us/call
    hk_conway_previous              0.628 ms   0.38 Mcyc   1.0%
      x0.44  1435.933 us/call
    is_build_draw                  38.923 ms  23.35 Mcyc  63.6%
      x1.00  38923.431 us/call
      is_build_scan                37.501 ms  22.50 Mcyc  61.3%
        x1.00  37500.732 us/call
        scan_face_setup             0.457 ms   0.27 Mcyc   0.7%
          x15.88  28.804 us/call
        mesh_copy_capture           0.239 ms   0.14 Mcyc   0.4%
          x18.25  13.078 us/call
        mesh_capture                0.829 ms   0.50 Mcyc   1.4%
          x3.62  228.697 us/call
        scan_face_sweep             1.768 ms   1.06 Mcyc   2.9%
          x1.44  1229.634 us/call
          scan_sweep_column         0.561 ms   0.34 Mcyc   0.9%
            x121.94  4.598 us/call
        mesh_sweep_paths            0.044 ms   0.03 Mcyc   0.1%
          x1.44  30.607 us/call
        mesh_classify               0.516 ms   0.31 Mcyc   0.8%
          x0.75  688.455 us/call
      is_mesh_transform             0.078 ms   0.05 Mcyc   0.1%
        x1.00  78.115 us/call
    hk_conway_compile               0.391 ms   0.23 Mcyc   0.6%
      x1.00  390.919 us/call
    hk_conway_sweep                 0.902 ms   0.54 Mcyc   1.5%
      x1.00  901.653 us/call
  is_ripple_prepare                 0.003 ms   0.00 Mcyc   0.0%
    x1.00  2.520 us/call
  canvas_clear                      0.094 ms   0.06 Mcyc   0.2%
    x1.00  94.100 us/call
  canvas_buffer_wait               18.946 ms  11.37 Mcyc  31.0%
    x1.00  18945.679 us/call
```

Wall min/avg/max = 39.932/61.221/73.096 ms. Render window mean is 42.261 ms; its peak frame is 50.202 ms. The difference from wall includes display-sync idle.

`is_timeline_step`: 1.00 calls/frame, 41778.812 µs/call. Scope tags, where present, retain the parser limitations.

### Peak-render window (window frames 241–256)

```text
scan_mesh_raster                   26.569 ms  15.94 Mcyc  42.1%
  x224.19  118.511 us/call
scan_face_setup                     2.432 ms   1.46 Mcyc   3.9%
  x237.19  10.254 us/call
mesh_visibility                     0.499 ms   0.30 Mcyc   0.8%
  x1.00  498.945 us/call
frame                              63.131 ms  37.88 Mcyc 100.0%
  x1.00  63131.179 us/call
  pov_preserve_half                 0.140 ms   0.08 Mcyc   0.2%
    x1.00  139.891 us/call
  is_timeline_step                 40.506 ms  24.30 Mcyc  64.2%
    x1.00  40506.313 us/call
    hk_conway_previous              0.566 ms   0.34 Mcyc   0.9%
      x0.56  1005.732 us/call
    is_build_draw                  37.685 ms  22.61 Mcyc  59.7%
      x1.00  37685.173 us/call
      is_build_scan                36.372 ms  21.82 Mcyc  57.6%
        x1.00  36371.794 us/call
        scan_face_setup             0.448 ms   0.27 Mcyc   0.7%
          x13.81  32.460 us/call
        mesh_copy_capture           0.209 ms   0.13 Mcyc   0.3%
          x14.75  14.181 us/call
        mesh_capture                0.739 ms   0.44 Mcyc   1.2%
          x3.56  207.499 us/call
        scan_face_sweep             2.126 ms   1.28 Mcyc   3.4%
          x0.69  3092.349 us/call
          scan_sweep_column         0.630 ms   0.38 Mcyc   1.0%
            x78.19  8.057 us/call
        mesh_sweep_paths            0.024 ms   0.01 Mcyc   0.0%
          x0.69  34.572 us/call
        mesh_classify               0.456 ms   0.27 Mcyc   0.7%
          x0.81  561.326 us/call
      is_mesh_transform             0.074 ms   0.04 Mcyc   0.1%
        x1.00  73.734 us/call
    hk_conway_compile               0.350 ms   0.21 Mcyc   0.6%
      x1.00  349.841 us/call
    hk_conway_sweep                 0.828 ms   0.50 Mcyc   1.3%
      x1.00  828.015 us/call
  is_ripple_prepare                 0.002 ms   0.00 Mcyc   0.0%
    x1.00  2.409 us/call
  canvas_clear                      0.095 ms   0.06 Mcyc   0.2%
    x1.00  94.729 us/call
  canvas_buffer_wait               22.144 ms  13.29 Mcyc  35.1%
    x1.00  22143.767 us/call
```

Wall min/avg/max = 48.624/63.145/73.120 ms. Render window mean is 40.987 ms; its peak frame is 53.854 ms. The difference from wall includes display-sync idle.

`is_timeline_step`: 1.00 calls/frame, 40506.313 µs/call. Scope tags, where present, retain the parser limitations.

### Per-preset table

The capture has 16 distinct shape markers and returns to its first shape, confirming a closed cycle. Ownership peaks include transitions and exclude setup frame 1. Clean windows use the parser single-owner/modal-call selection; continuously constructed shapes may lack a static hold.

| Shape | V/E/F/I at spawn | Blends/frame | Clean timeline ms | Clean render ms | Clean fps | Peak render ms | Baseline peak ms | Δ peak ms | Spilled/frames | Clean/all windows |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| truncatedCube_hk90_meta | 24/36/14/72 | unavailable | 41.779 | 42.261 | 16.00 | 53.854 | 53.491 | +0.363 | 0/276 | 14/17 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 20/30/12/60 | unavailable | 40.978 | 41.290 | 16.00 | 53.665 | 53.634 | +0.031 | 0/264 | 13/16 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 120/180/62/360 | unavailable | 38.162 | 39.150 | 16.00 | 51.447 | 51.484 | -0.037 | 0/234 | 10/13 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 60/90/32/180 | unavailable | 35.437 | 35.752 | 16.00 | 45.397 | 44.929 | +0.468 | 0/216 | 10/13 |
| snubDodecahedron_truncate5d_ambo_dual | 60/150/92/300 | unavailable | 32.166 | 33.260 | 16.00 | 36.397 | 36.185 | +0.212 | 0/234 | 12/15 |
| dodecahedron_bevel2_relax_gyro | 20/30/12/60 | unavailable | 30.282 | 30.832 | 16.00 | 35.040 | 34.967 | +0.073 | 0/264 | 13/16 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 12/30/20/60 | unavailable | 29.004 | 29.319 | 16.00 | 42.300 | 41.733 | +0.567 | 0/160 | 8/10 |
| icosahedron_snub_relax_truncate033_hankin62 | 12/30/20/60 | unavailable | 27.988 | 28.297 | 16.00 | 38.543 | 38.236 | +0.307 | 0/144 | 7/9 |
| rhombicuboctahedron_hk63_ambo_hk63 | 24/48/26/96 | unavailable | 27.480 | 27.784 | 16.00 | 36.596 | 36.642 | -0.046 | 0/177 | 9/11 |
| dodecahedron_hk72_ambo_dual_hk20 | 20/30/12/60 | unavailable | 27.115 | 27.420 | 16.00 | 34.718 | 34.805 | -0.087 | 0/176 | 9/11 |
| truncatedIcosahedron_hk58_chamfer63 | 60/90/32/180 | unavailable | 26.772 | 27.084 | 16.00 | 34.750 | 34.769 | -0.019 | 0/186 | 9/12 |
| icosidodecahedron_truncate5d_ambo_dual | 30/60/32/120 | unavailable | 26.325 | 26.731 | 16.00 | 35.135 | 34.594 | +0.541 | 0/234 | 12/15 |
| dodecahedron_hk54_ambo_hk72 | 20/30/12/60 | unavailable | 24.455 | 24.762 | 16.00 | 30.965 | 30.981 | -0.016 | 0/140 | 6/8 |
| dodecahedron_hk62_ambo_hk62 | 20/30/12/60 | unavailable | 23.707 | 24.015 | 16.00 | 31.145 | 31.332 | -0.187 | 0/208 | 11/14 |
| icosahedron_ambo_truncate033_hankin59 | 12/30/20/60 | unavailable | 23.575 | 23.886 | 16.00 | 35.938 | 35.802 | +0.136 | 0/204 | 11/14 |
| octahedron_hk17_ambo_hk73 | 6/12/8/24 | unavailable | 22.835 | 23.144 | 16.00 | 37.088 | 36.804 | +0.284 | 0/210 | 11/14 |

### Per-pixel figures

Blended-pixel counters are disabled in this timing image. Blend coverage, cycles/blend and scan cycles/blend are unavailable; geometric quadrant size is not a measured blended-pixel count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; ISR interval 1012.586 ms.

```text
isr_wake       1166.56/frame
  min/avg/max 0.541/2.005/27.165 us  cpu 3.69%
isr_pack        145.81/frame
  min/avg/max 8.741/9.423/11.896 us  cpu 2.17%
isr_dma_submit  145.81/frame
  min/avg/max 0.670/0.935/9.068 us  cpu 0.21%
```

- Submission is CPU setup time; packing performs the larger measured CPU work. Asynchronous LED wire-transfer duration is not measured.
- Inclusive `isr_wake` consumes 3.69% of its ISR interval, leaving 60.194 ms per display window outside that measured interrupt work. Packing and submission are nested in wake; their shares must not be added.
- Render scopes already include interrupts. DMA-completion and other interrupt costs are not separately measured.

## Summary ranking

1. `is_timeline_step` — 64.16% of this root window, 40.506 ms/frame.
2. `is_build_draw` — 59.69% of this root window, 37.685 ms/frame.
3. `is_build_scan` — 57.61% of this root window, 36.372 ms/frame.

Matched COM4 / TS4 comparison; both use the same 210-second capture budget and shipping selective-O3 configuration.

| Capture | Steady peak ms | Steady mean ms | Spilled/frames |
|---|--:|--:|--:|
| Before `d459cc1969e4` | 53.634 | 25.009 | 0/3327 |
| After `771290aa484e` | 53.854 | 25.159 | 0/3327 |
| Δ after − before | +0.220 | +0.150 | +0 spilled |

Matched schedule comparison uses common frame numbers 2–3328 with identical per-frame shape ownership (3327 frames; 0 ownership mismatches excluded). Frame 1 is excluded from both. This removes capture-tail weighting differences.

| Common schedule frames | Peak ms | Mean ms | Spilled/frames |
|---|--:|--:|--:|
| Before | 53.634 | 25.009 | 0/3327 |
| After | 53.854 | 25.159 | 0/3327 |
| Δ after − before | +0.220 | +0.150 | +0 spilled |

| Shape, matched schedule frames | Before mean ms | After mean ms | Δ mean ms | Frames |
|---|--:|--:|--:|--:|
| dodecahedron_bevel2_relax_gyro | 21.504 | 21.633 | +0.129 | 264 |
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 32.199 | 32.344 | +0.144 | 264 |
| dodecahedron_hk54_ambo_hk72 | 22.298 | 22.418 | +0.120 | 140 |
| dodecahedron_hk62_ambo_hk62 | 21.455 | 21.561 | +0.106 | 208 |
| dodecahedron_hk72_ambo_dual_hk20 | 23.238 | 23.361 | +0.123 | 176 |
| icosahedron_ambo_truncate033_hankin59 | 20.444 | 20.555 | +0.111 | 204 |
| icosahedron_snub_relax_truncate033_hankin62 | 24.919 | 25.021 | +0.102 | 144 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 26.270 | 26.394 | +0.124 | 160 |
| icosidodecahedron_truncate5d_ambo_dual | 19.515 | 19.698 | +0.182 | 234 |
| octahedron_hk17_ambo_hk73 | 20.375 | 20.510 | +0.135 | 210 |
| rhombicuboctahedron_hk63_ambo_hk63 | 24.212 | 24.329 | +0.117 | 177 |
| snubDodecahedron_truncate5d_ambo_dual | 23.483 | 23.626 | +0.144 | 234 |
| truncatedCube_hk90_meta | 31.142 | 31.398 | +0.256 | 276 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 30.019 | 30.228 | +0.209 | 216 |
| truncatedIcosahedron_hk58_chamfer63 | 24.531 | 24.677 | +0.147 | 186 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 30.437 | 30.608 | +0.171 | 234 |

Systematic root-versus-wall diagnostic: baseline 177.5–259.9 ppm across 208 windows; 0 windows meet 5 ppm; final 193.8–270.2 ppm across 208 windows; 0 windows meet 5 ppm. Positive differences mean wall exceeds inner root cycles. These diagnostics remain separate from frame-telemetry cadence results.

| Image | RAM1 code before/after/Δ B | RAM1 variables before/after/Δ B | FLASH data before/after/Δ B |
|---|--:|--:|--:|
| profile | 30,792/30,840/+48 | 313,664/313,664/+0 | 208,804/208,884/+80 |
| phantasm | 176,024/176,040/+16 | 312,736/312,736/+0 | 853,152/853,320/+168 |

README cells: peak 🟢 53.854 (16), spilled 🟢 0/3327 (0.00%).

Global-O3 captures retain their original source and are unchanged.

## Caveats

- All cycle scopes absorb ISR time and instrumentation overhead; instrumenting per-pixel work can change its cost.
- `filter_blend` parenting can hide its subtree when its first parent has no calls; no measured blend count is claimed here.
- The shipping image uses selective `HS_O3` transforms and scan kernels, not a global-O3 image.
- TS4 divides every stage duration, changing the sampled animation workload; both captures use TS4. `HS_PROFILE_EPOCH_REVS=1920` shortens the epoch to 240 s, still longer than the 210 s capture.
- Both source snapshots were committed worker/integrator trees; capture provenance and adjacent build/envdump files identify their images. Local ignored build outputs do not imply source changes.
- Frame telemetry supports cadence comparison; strict cycle/wall exactness remains subject to the diagnostics below. Single paired captures do not establish statistical significance for small timing deltas. The observed +0.220 ms peak and +0.150 ms mean are increases in this pair; they are not evidence of a speedup.

Baseline parser validation:

```text
effect=IslamicStars  windows=208
  [PASS] captured >=3 windows (208)
  [PASS] every window names one effect and resolution (IslamicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] shape markers advance and return to the first (cycle closed): 16 distinct
  [PASS] per-frame render is render, not wall (0 of 208 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 208 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): is_mesh_transform, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, mesh_visibility, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'is_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (199.3 ppm, frames 1473-1488)

INVALID: IslamicStars
```

Final parser validation:

```text
effect=IslamicStars  windows=208
  [PASS] captured >=3 windows (208)
  [PASS] every window names one effect and resolution (IslamicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] shape markers advance and return to the first (cycle closed): 16 distinct
  [PASS] per-frame render is render, not wall (0 of 208 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 208 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): is_mesh_transform, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, mesh_visibility, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'is_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (239.4 ppm, frames 1473-1488)

INVALID: IslamicStars
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`; the [Setup command](#setup) builds, flashes and captures this full-cycle TS4 comparison through the board lock.
