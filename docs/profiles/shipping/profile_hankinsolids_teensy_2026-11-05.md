# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HankinSolids`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/hankinsolids_ship.log`. Captured 2026-11-04 23:22 America/Los_Angeles on COM4. Replaces profile_hankinsolids_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 170 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh HankinSolids profile 170 16` |

Image size:

```text
teensy_size:   FLASH: code:160720, data:164012, headers:9092   free for files:1697792
teensy_size:    RAM1: variables:314656, code:31448, padding:1320   free for local variables:176864
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `2609ca19ede9a7d1816d7c2aa286803aae88608662411bbc032d9b486cd1f0cb`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1409–1424, root cyc/600 versus measured wall sum differs by **0.413 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   HankinSolids             worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 38,564 B  full-set worst 4
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **30.408 ms** (frame 2385), **0/2687 spilled** (0.00%). Setup frame 1: 19.120 ms. Complete windows cover frames 1–2688; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 16.901 ms/frame. `hk_timeline_step` averages 16.630 ms/frame over complete windows; its worst window is 26.924 ms/frame (frames 1409–1424).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 8.881–30.408 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     56.736 ms 34.042 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.082 Mcyc   0.2%
    x0.94/frame; 144.998 us/call
  hk_timeline_step                        15.916 ms  9.550 Mcyc  28.1%
    hk_update_hankin                       0.071 ms  0.043 Mcyc   0.1%
      x0.94/frame; 75.727 us/call
    hk_draw_mesh                          15.628 ms  9.377 Mcyc  27.5%
      hk_mesh_scan                        15.338 ms  9.203 Mcyc  27.0%
        scan_face_sweep                    7.751 ms  4.650 Mcyc  13.7%
          scan_sweep_column                2.015 ms  1.209 Mcyc   3.6%
            x594.00/frame; 3.393 us/call
        mesh_sweep_paths                   0.033 ms  0.020 Mcyc   0.1%
          x4.12/frame; 7.888 us/call
        scan_face_setup                    0.143 ms  0.086 Mcyc   0.3%
          x4.00/frame; 35.653 us/call
        mesh_copy_capture                  0.036 ms  0.022 Mcyc   0.1%
          x4.00/frame; 9.010 us/call
        mesh_capture                       0.009 ms  0.005 Mcyc   0.0%
          x1.88/frame; 4.728 us/call
        mesh_classify                      0.092 ms  0.055 Mcyc   0.2%
          x0.81/frame; 113.167 us/call
        scan_mesh_raster                   6.750 ms  4.050 Mcyc  11.9%
          x6.25/frame; 1080.058 us/call
        scan_face_setup                    0.064 ms  0.038 Mcyc   0.1%
          x2.25/frame; 28.338 us/call
      hk_mesh_transform                    0.104 ms  0.063 Mcyc   0.2%
        x1.00/frame; 104.177 us/call
    hk_update_hankin                       0.073 ms  0.044 Mcyc   0.1%
      x1.00/frame; 72.931 us/call
  canvas_clear                             0.101 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.595 us/call
  canvas_buffer_wait                      40.547 ms 24.328 Mcyc  71.5%
    x1.00/frame; 40546.965 us/call
```

Wall min/avg/max: 16.591/56.735/64.565 ms. Window mean render is 16.189 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 2385–2400)

```text
frame                                     61.765 ms 37.059 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 144.687 us/call
  hk_timeline_step                        24.542 ms 14.725 Mcyc  39.7%
    hk_update_hankin                       0.671 ms  0.403 Mcyc   1.1%
      x1.00/frame; 670.905 us/call
    hk_draw_mesh                          23.046 ms 13.828 Mcyc  37.3%
      hk_mesh_scan                        22.248 ms 13.349 Mcyc  36.0%
        scan_face_setup                    0.123 ms  0.074 Mcyc   0.2%
          x5.88/frame; 20.902 us/call
        mesh_copy_capture                  0.021 ms  0.012 Mcyc   0.0%
          x6.00/frame; 3.437 us/call
        mesh_capture                       0.190 ms  0.114 Mcyc   0.3%
          x5.00/frame; 38.069 us/call
        mesh_classify                      0.237 ms  0.142 Mcyc   0.4%
          x1.00/frame; 237.360 us/call
        scan_mesh_raster                  19.203 ms 11.522 Mcyc  31.1%
          x93.12/frame; 206.212 us/call
        scan_face_setup                    1.353 ms  0.812 Mcyc   2.2%
          x94.31/frame; 14.341 us/call
      hk_mesh_transform                    0.570 ms  0.342 Mcyc   0.9%
        x1.00/frame; 570.104 us/call
    hk_update_hankin                       0.693 ms  0.416 Mcyc   1.1%
      x1.00/frame; 693.216 us/call
  canvas_clear                             0.093 ms  0.056 Mcyc   0.1%
    x1.00/frame; 92.568 us/call
  canvas_buffer_wait                      36.950 ms 22.170 Mcyc  59.8%
    x1.00/frame; 36950.339 us/call
```

Wall min/avg/max: 60.366/61.765/62.468 ms. Window mean render is 24.815 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 1409–1424)

```text
frame                                     62.276 ms 37.366 Mcyc 100.0%
  pov_preserve_half                        0.146 ms  0.087 Mcyc   0.2%
    x1.00/frame; 145.637 us/call
  hk_timeline_step                        26.924 ms 16.154 Mcyc  43.2%
    hk_update_hankin                       0.302 ms  0.181 Mcyc   0.5%
      x1.00/frame; 301.556 us/call
    hk_draw_mesh                          26.153 ms 15.692 Mcyc  42.0%
      hk_mesh_scan                        25.702 ms 15.421 Mcyc  41.3%
        scan_face_sweep                    1.942 ms  1.165 Mcyc   3.1%
          scan_sweep_column                0.506 ms  0.303 Mcyc   0.8%
            x72.00/frame; 7.024 us/call
        mesh_sweep_paths                   0.009 ms  0.006 Mcyc   0.0%
          x0.50/frame; 18.697 us/call
        scan_face_setup                    0.062 ms  0.037 Mcyc   0.1%
          x1.88/frame; 33.043 us/call
        mesh_copy_capture                  0.014 ms  0.008 Mcyc   0.0%
          x1.88/frame; 7.404 us/call
        mesh_capture                       0.132 ms  0.079 Mcyc   0.2%
          x3.75/frame; 35.174 us/call
        mesh_classify                      0.137 ms  0.082 Mcyc   0.2%
          x0.75/frame; 182.361 us/call
        scan_mesh_raster                  21.967 ms 13.180 Mcyc  35.3%
          x47.56/frame; 461.854 us/call
        scan_face_setup                    0.751 ms  0.451 Mcyc   1.2%
          x46.38/frame; 16.199 us/call
      hk_mesh_transform                    0.256 ms  0.154 Mcyc   0.4%
        x1.00/frame; 255.990 us/call
    hk_update_hankin                       0.320 ms  0.192 Mcyc   0.5%
      x1.00/frame; 320.096 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.108 us/call
  canvas_buffer_wait                      35.081 ms 21.048 Mcyc  56.3%
    x1.00/frame; 35080.521 us/call
```

Wall min/avg/max: 55.869/62.276/68.809 ms. Window mean render is 27.196 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 2513–2528)

```text
frame                                     63.027 ms 37.816 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 145.176 us/call
  hk_timeline_step                        20.683 ms 12.410 Mcyc  32.8%
    hk_conway_previous                     0.042 ms  0.025 Mcyc   0.1%
      x0.31/frame; 134.809 us/call
    hk_draw_mesh                           3.454 ms  2.072 Mcyc   5.5%
      hk_mesh_scan                         3.397 ms  2.038 Mcyc   5.4%
        x0.31/frame; 10869.144 us/call
      hk_mesh_transform                    0.005 ms  0.003 Mcyc   0.0%
        x0.31/frame; 14.771 us/call
    hk_conway_compile                      0.023 ms  0.014 Mcyc   0.0%
      x0.31/frame; 72.075 us/call
    hk_conway_sweep                        0.050 ms  0.030 Mcyc   0.1%
      x0.31/frame; 160.009 us/call
    hk_update_hankin                       0.324 ms  0.194 Mcyc   0.5%
      x0.62/frame; 518.558 us/call
    hk_draw_mesh                          15.955 ms  9.573 Mcyc  25.3%
      hk_mesh_scan                        15.805 ms  9.483 Mcyc  25.1%
        scan_mesh_raster                  17.886 ms 10.731 Mcyc  28.4%
          x57.25/frame; 312.412 us/call
        scan_face_setup                    0.802 ms  0.481 Mcyc   1.3%
          x57.25/frame; 14.004 us/call
      hk_mesh_transform                    0.019 ms  0.011 Mcyc   0.0%
        x0.69/frame; 27.419 us/call
    hk_update_hankin                       0.366 ms  0.219 Mcyc   0.6%
      x0.69/frame; 531.966 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 86.286 us/call
  canvas_buffer_wait                      42.088 ms 25.253 Mcyc  66.8%
    x1.00/frame; 42088.322 us/call
```

Wall min/avg/max: 55.869/63.027/75.191 ms. Window mean render is 20.939 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| truncatedCuboctahedron truncatedCuboctahedron | 8/7 | unavailable | 26.924 | 27.196 | 16.06 | modal-call-count hold |
| truncatedIcosidodecahedron truncatedIcosidodecahedron | 7/6 | unavailable | 26.589 | 26.857 | 15.93 | modal-call-count hold |
| snubDodecahedron snubDodecahedron | 7/6 | unavailable | 26.123 | 26.388 | 15.94 | modal-call-count hold |
| truncatedDodecahedron truncatedDodecahedron | 7/6 | unavailable | 25.892 | 26.170 | 16.07 | modal-call-count hold |
| truncatedCube truncatedCube | 7/6 | unavailable | 25.617 | 25.892 | 16.02 | modal-call-count hold |
| rhombicosidodecahedron rhombicosidodecahedron | 7/6 | unavailable | 25.479 | 25.744 | 15.96 | modal-call-count hold |
| truncatedIcosahedron truncatedIcosahedron | 7/6 | unavailable | 24.790 | 25.056 | 16.05 | modal-call-count hold |
| truncatedOctahedron truncatedOctahedron | 7/6 | unavailable | 22.485 | 22.759 | 15.99 | modal-call-count hold |
| snubCube snubCube | 8/7 | unavailable | 22.416 | 22.684 | 15.99 | modal-call-count hold |
| icosidodecahedron icosidodecahedron | 17/10 | unavailable | 21.607 | 21.872 | 15.98 | modal-call-count hold |
| rhombicuboctahedron rhombicuboctahedron | 7/6 | unavailable | 20.316 | 20.585 | 16.08 | modal-call-count hold |
| truncatedTetrahedron truncatedTetrahedron | 6/5 | unavailable | 20.176 | 20.454 | 15.99 | modal-call-count hold |
| dodecahedron dodecahedron | 15/10 | unavailable | 19.091 | 19.358 | 15.96 | modal-call-count hold |
| cuboctahedron cuboctahedron | 16/9 | unavailable | 18.462 | 18.734 | 16.03 | modal-call-count hold |
| cube cube | 13/8 | unavailable | 17.954 | 18.227 | 16.13 | modal-call-count hold |
| octahedron octahedron | 11/8 | unavailable | 17.937 | 18.222 | 16.06 | modal-call-count hold |
| icosahedron icosahedron | 9/6 | unavailable | 17.814 | 18.085 | 15.97 | modal-call-count hold |
| 0  | 7/6 | unavailable | 15.916 | 16.189 | n/a | startup window (includes setup frame 1) |
| tetrahedron tetrahedron | 2/1 | unavailable | 12.293 | 12.575 | 16.00 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| truncatedIcosidodecahedron | 30.408 | 0/117 | 16.00 |
| truncatedDodecahedron | 30.008 | 0/105 | 16.00 |
| snubDodecahedron | 29.795 | 0/117 | 16.00 |
| truncatedCuboctahedron | 29.268 | 0/117 | 16.00 |
| truncatedCube | 28.799 | 0/105 | 16.00 |
| snubCube | 28.360 | 0/117 | 16.00 |
| rhombicosidodecahedron | 28.092 | 0/117 | 16.00 |
| truncatedIcosahedron | 26.854 | 0/105 | 16.00 |
| truncatedTetrahedron | 25.162 | 0/105 | 16.00 |
| icosidodecahedron | 25.132 | 0/285 | 16.00 |
| truncatedOctahedron | 24.533 | 0/105 | 16.00 |
| rhombicuboctahedron | 23.391 | 0/105 | 16.00 |
| cube | 22.110 | 0/213 | 16.00 |
| icosahedron | 21.911 | 0/160 | 16.00 |
| dodecahedron | 21.891 | 0/225 | 16.00 |
| octahedron | 21.185 | 0/177 | 16.00 |
| cuboctahedron | 21.116 | 0/285 | 16.00 |
| 0 | 17.978 | 0/103 | 16.00 |
| tetrahedron | 15.243 | 0/24 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.00/frame
  min/avg/max 0.553/1.963/25.436 us; cpu 3.61%
  isr_pack         143.88/frame
    min/avg/max 8.806/9.440/11.943 us; cpu 2.17%
  isr_dma_submit   143.88/frame
    min/avg/max 0.691/0.928/0.998 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999056 us
isr_wake       n=18416 cyc min/avg/max=332/1178/15262 ns min/avg/max=553/1963/25436 total=36159 us cpu=3.61%
isr_pack       n=2302 cyc min/avg/max=5284/5664/7166 ns min/avg/max=8806/9440/11943 total=21733 us cpu=2.17%
isr_dma_submit n=2302 cyc min/avg/max=415/557/599 ns min/avg/max=691/928/998 total=2140 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.61% leaves approximately 60.244 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `hk_timeline_step` — 24.542 ms/frame, 39.73% of this root window.
2. `hk_draw_mesh` — 23.046 ms/frame, 37.31% of this root window.
3. `hk_mesh_scan` — 22.248 ms/frame, 36.02% of this root window.
4. `scan_mesh_raster` — 19.203 ms/frame, 31.09% of this root window.
5. `scan_face_setup` — 1.353 ms/frame, 2.19% of this root window.

README cells: peak 🟢 30.408 (19), spilled 🟢 0/2687 (0.00%).

README setup frame 1: 19.120 ms. Captured 2026-11-04 23:22 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=HankinSolids  windows=168
  [PASS] captured >=3 windows (168)
  [PASS] every window names one effect and resolution (HankinSolids 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] shape markers advance and return to the first (cycle closed): 18 distinct
  [PASS] per-frame render is render, not wall (0 of 168 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 168 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): hk_draw_mesh, hk_mesh_scan, hk_mesh_transform, hk_update_hankin, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'hk_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.4 ppm, frames 1409-1424)

VALID: HankinSolids
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HankinSolids`, `HS_PROFILE_WINDOW=16`. `just profile HankinSolids` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
