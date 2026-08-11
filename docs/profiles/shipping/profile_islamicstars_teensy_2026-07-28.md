# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-07-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `build/prof/islamicstars_ship.log`. Replaces
`profile_islamicstars_teensy_2026-07-27.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, board COM4 |
| Image | `profile` env: `-Os` + newlib-nano + DMA LEDs + `HS_PROFILE_ENABLE`; the per-face SDF raster driver (R1/R2), `SDF::Face::plane_dsq_sector`, `ripple_transform` and the effect's own `transform_shape`/`draw_build_mesh` run inside `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `542a5b49` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, `HS_PROFILE_EPOCH_REVS=1920`, `HS_PROFILE_TRANS_SPEED=4` |
| Reproduce | `just profile IslamicStars` |

Image size: `FLASH: code:126060, data:577368, headers:8244` / `RAM1:
variables:314848, code:44104, padding:21432, free:143904` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 2993–3008 root counter cyc ÷ 600 MHz
(624,566,407 ÷ 600 = 1,040,944 µs) matches the measured wall sum
(1,040,944 µs) within **0.0 ppm** (`tools/parse_profile.py ... validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `is_timeline_step`
18.7–47.1 ms/f across the 24 shapes, worst hold
`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` at 47.05 ms/f, peak frame render
**56.91 ms** (frames 2993–3008), spilled **0/3328 frames (0%)**.

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px. Every
one of the 24 shapes holds 16 fps, and so does every build leg — the peak frame
of the pass is a build window (the `truncatedIcosidodecahedron_truncate50d_ambo_dual`
dual-bridge legs), not a hold window, and it still lands inside one display
window with 5.6 ms to spare. The `canvas_buffer_wait` scope is the round-up idle
to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the effect alternates between holding a shape (steady per-face
SDF raster) and an opchain build leg that compiles the next Conway operator
chain while still drawing. 43 `Spawning Shape:` markers over the pass, all 24
shapes visited, cycle closed.

### Shape hold — worst of the capture (window frames 897–912)

```
frame                   62.50 ms  37.50 Mcyc  100%
  is_timeline_step      47.05 ms  28.23 Mcyc   75%
    is_draw_shape       47.00 ms  28.20 Mcyc   75%  x1
      is_mesh_scan      43.00 ms  25.80 Mcyc   69%  x1
        scan_mesh_ras.  32.63 ms  19.58 Mcyc   52%  x1082  30.2 us/call
          filter_blend   1.18 ms   0.71 Mcyc    2%  x18165  39.0 cyc/blend
        scan_face_setup 10.15 ms   6.09 Mcyc   16%  x1082   9.4 us/call
      is_mesh_transform  3.50 ms   2.10 Mcyc    6%  x1
      is_face_offsets    0.49 ms   0.30 Mcyc    1%  x1
  is_ripple_prepare      6.69 us   4048 cyc     0%
  canvas_clear          88.13 us  52880 cyc     0%
  canvas_buffer_wait    15.36 ms   9.21 Mcyc   25%
```

Wall min/avg/max = 61.74/62.50/63.39 ms — the tightest wall spread in the
capture, because a hold does identical work every frame. This is the costliest
hold on the roster: the `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` chain
builds out to 1,082 faces, so `scan_mesh_raster` is entered 1,082×/frame and
`scan_face_setup`'s per-face cost (9.4 µs) alone reaches 10.2 ms. Render is
47.1 ms against the 62.5 ms window, leaving 25% of the frame as display-sync
idle — the narrowest hold margin of the pass, and still one window.

### Opchain build leg — peak frame of the capture (window frames 2993–3008)

```
frame                   65.06 ms  39.04 Mcyc  100%
  is_timeline_step      30.78 ms  18.47 Mcyc   47%
    is_build_draw       25.88 ms  15.53 Mcyc   40%  x0.75  34.5 ms/call
      is_build_scan     25.85 ms  15.51 Mcyc   40%  x0.75
      is_mesh_transform 31.38 us  18842 cyc     0%  x0.75
    hk_conway_compile   97.38 us  58428 cyc     0%  x0.75
    hk_conway_op         0.53 ms   0.32 Mcyc    1%  x0.75
    is_draw_shape        2.78 ms   1.67 Mcyc    4%  x0.25
      is_mesh_scan       2.77 ms   1.66 Mcyc    4%  x0.25
        scan_mesh_ras.  26.22 ms  15.73 Mcyc   40%  x287   91.4 us/call
          filter_blend   0.98 ms   0.59 Mcyc    2%  x15724  37.2 cyc/blend
        scan_face_setup  2.35 ms   1.41 Mcyc    4%  x287    8.2 us/call
      is_face_offsets    6.44 us   3880 cyc     0%  x0.25
      is_mesh_transform  1.94 us   1194 cyc     0%  x0.25
  is_ripple_prepare      0.06 us     42 cyc     0%
  canvas_clear          86.63 us  52003 cyc     0%
  canvas_buffer_wait    34.19 ms  20.51 Mcyc   53%
```

Wall min/avg/max = 43.98/65.06/76.83 ms — this spread is **not** cadence jitter;
see the wall-metric caveat below. Build frames and draw frames alternate inside
the window (12 build calls against 4 draw calls across the 16 frames), and wall
tracks the frame-to-frame *change* in render, so alternating cost produces an
alternating wall while every frame still occupies exactly one display window.

`scan_mesh_raster`
reads 40% of the frame against `is_mesh_scan`'s 4%: it first entered under
`is_mesh_scan` and parents there, but the build path enters it too and the
counter carries both — its 287 calls/frame and 91.4 µs/call are the real
figures, the percentage against its parent is the artifact. `hk_conway_compile`
is 0.1 ms and `hk_conway_op` 0.5 ms: compiling and evaluating the operator
chain are free next to rasterizing alongside them.

### Per-preset table

Read from each shape's clean-hold windows (`parse_profile.py ... presets`);
build-straddle windows are excluded by the clean-hold criterion. All 24 shapes
visited and the cycle closed (24 distinct shapes over 43 spawn markers).
Windows per shape: 4–14. `V/E/F/I` is the **seed** solid from the
`Spawning Shape:` line — the operator chain expands it substantially before it
is drawn (the top row's seed dodecahedron builds out to V=3240, F=1082).
Ranked by dominant-scope cost.

The last column is that shape's peak *frame* render including the transition
that follows it (`parse_profile.py ... buckets`), which is the stricter figure
the aggregate README ranks on; the `is_timeline_step` column is the clean hold.

| # | shape | V/E/F/I | is_timeline_step ms | peak render ms | fps |
|---|---|---|--:|--:|--:|
| 1 | dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 20/30/12/60 | 47.05 | 49.51 | 16 |
| 2 | truncatedIcosahedron_ambo_relax_truncate001_hankin73 | 60/90/32/180 | 41.10 | 46.57 | 16 |
| 3 | truncatedIcosidodecahedron_bevel5_relax_hk77 | 120/180/62/360 | 40.50 | 46.72 | 16 |
| 4 | truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 60/90/32/180 | 39.11 | 41.27 | 16 |
| 5 | truncatedIcosahedron_ambo_relax100_hk54_needle | 60/90/32/180 | 36.06 | 53.88 | 16 |
| 6 | truncatedOctahedron_gyro_kis_hk17 | 24/36/14/72 | 34.93 | 45.35 | 16 |
| 7 | truncatedIcosidodecahedron_truncate50d_ambo_dual | 120/180/62/360 | 34.91 | 56.91 | 16 |
| 8 | truncatedIcosahedron_hk54_ambo_hk72 | 60/90/32/180 | 33.96 | 38.71 | 16 |
| 9 | truncatedIcosahedron_ambo_relax_truncate33_hk64 | 60/90/32/180 | 32.33 | 37.27 | 16 |
| 10 | truncatedIcosahedron_truncate50d_ambo_dual | 60/90/32/180 | 30.09 | 44.39 | 16 |
| 11 | dodecahedron_ambo_bevel33_relax_hk66 | 20/30/12/60 | 29.51 | 32.94 | 16 |
| 12 | icosahedron_snub_relax_truncate033_hankin62 | 12/30/20/60 | 29.33 | 30.03 | 16 |
| 13 | truncatedIcosahedron_hk58_chamfer63 | 60/90/32/180 | 28.34 | 29.32 | 16 |
| 14 | snubDodecahedron_truncate5d_ambo_dual | 60/150/92/300 | 28.14 | 37.95 | 16 |
| 15 | dodecahedron_bevel2_relax_gyro | 20/30/12/60 | 27.73 | 37.89 | 16 |
| 16 | rhombicuboctahedron_hk63_ambo_hk63 | 24/48/26/96 | 25.07 | 28.42 | 16 |
| 17 | dodecahedron_hk72_ambo_dual_hk20 | 20/30/12/60 | 23.52 | 30.06 | 16 |
| 18 | icosahedron_kis_gyro | 12/30/20/60 | 23.13 | 30.67 | 16 |
| 19 | dodecahedron_hk54_ambo_hk72 | 20/30/12/60 | 22.93 | 24.82 | 16 |
| 20 | dodecahedron_hk62_ambo_hk62 | 20/30/12/60 | 21.76 | 24.80 | 16 |
| 21 | icosahedron_ambo_truncate033_hankin59 | 12/30/20/60 | 21.57 | 24.03 | 16 |
| 22 | octahedron_hk17_ambo_hk73 | 6/12/8/24 | 20.56 | 24.08 | 16 |
| 23 | icosidodecahedron_truncate5d_ambo_dual | 30/60/32/120 | 19.56 | 23.60 | 16 |
| 24 | octahedron_hk34_ambo_hk72 | 6/12/8/24 | 18.74 | 21.44 | 16 |

The span is 2.5× (18.7–47.1 ms) and does not track seed size: the 20-vertex
dodecahedron chain tops the table while a 120-vertex truncated
icosidodecahedron sits mid-pack. Post-operator face count and the resulting
probe density set the cost, not the seed mesh.

### Per-pixel figures

Worst hold: 18,165 blended px/frame against the 10,368 px quadrant = **1.75×
coverage** (the star interiors overdraw), at **39.0 cyc/blend**.
`scan_mesh_raster` spends 32.63 ms/frame over those blends = 1,078 cyc per
blended pixel, so ~96% of the raster is probe evaluation that never reaches a
blend — consistent with the effect being probe-bound rather than blend-bound.
The build leg's ratio is the same shape (15,724 blends at 37.2 cyc/blend
against a 26.22 ms raster).

## Column-ISR / DMA marshaling cost

```
isr_wake        x1200/frame  min/avg/max 0.58/1.79/16.27 us  cpu 3.30%
isr_pack        x150/frame   min/avg/max 6.11/6.92/9.10 us   cpu 1.59%
isr_dma_submit  x150/frame   min/avg/max 0.62/0.93/1.06 us   cpu 0.21%
```

- Submit is 13% of pack's CPU cost — the column marshal dominates, the DMA kick
  is nearly free.
- The wire transfer is asynchronous; submit returns in 0.93 µs and the engine
  runs against the render.
- Total ISR share is 5.10%, leaving ≈ 59.3 ms of render budget per 62.5 ms
  window. The worst hold uses 47.1 ms and the peak build frame 56.9 ms, so both
  clear the budget without a speedup — the peak leaves 2.4 ms of margin.

## Summary ranking

1. `scan_mesh_raster` — 52% of the worst hold frame, 32.63 ms: per-face SDF
   probing over 1,082 faces, 1,078 cyc per blended pixel.
2. `is_build_scan` — 40% of the peak build frame, 25.85 ms: the opchain build
   leg's own raster, at 34.5 ms per build call the single most expensive scope
   in the capture.
3. `scan_face_setup` — 16% of the worst hold, 10.15 ms: per-face setup at
   9.4 µs/call over 1,082 faces/frame; its share grows with face count, so it
   is the second lever on the heavy chains.
4. `is_mesh_transform` — 6%, 3.50 ms: whole-mesh transform once per frame,
   only material on the high-face-count holds.
5. `hk_conway_op` + `hk_conway_compile` — 1% combined, 0.63 ms: operator
   evaluation and chain compilation are not the build leg's cost; the raster
   that runs alongside them is.

Against the ledger: the effect is probe-bound (~464 cyc/probe measured
2026-07-16), and the sector-walk raster plus the ripple amplitude fold are the
landed wins that put every shape inside one display window.

## Caveats

- **Wall min/avg/max is not a cadence signal.** `buffer_wait` sits at the head of
  `draw_frame`, so a measured frame runs from the end of the previous render to
  the end of this one and wall picks up the render *derivative*:
  `wall(N) = P + render(N) - render(N-1)`, which reproduces every frame of this
  capture to a 7 us median / 0.28 ms max residual (P = 62.46 ms). The display
  flip itself is rigid and `buffer_wait` absorbs the remainder exactly; total
  pass wall is 0.999x `frames x 62.5 ms`, i.e. one display window per frame with
  no doubled frames. A volatile per-frame cost therefore shows a wide wall
  spread while sitting perfectly on cadence — read `spilled` and peak render,
  never wall max.
- All scopes absorb ISR time (CYCCNT free-runs) — the 5.10% ISR share is inside
  every number above.
- `filter_blend` parents under `scan_mesh_raster`; in the build window its
  parent chain reports >100% of `is_mesh_scan` because `scan_mesh_raster` is
  entered from both the draw and build paths while parenting only to the first.
  Its calls ≈ blended pixels regardless.
- The per-face SDF raster driver (`Scan::Mesh::draw` R1/R2),
  `SDF::Face::plane_dsq_sector`, `ripple_transform`, and the effect's own
  `transform_shape` / `draw_build_mesh` sit inside `HS_O3` regions, which
  activate only on the `-Os` device build — the `profile` env therefore
  measures the shipping selective-O3 configuration.
- `HS_PROFILE_TRANS_SPEED=4` compresses carousel dwell so all 24 shapes fit one
  capture. It changes how long a shape *holds*, never its per-frame render cost,
  so the per-shape ms above are the shipping figures.
- Captured at tip `542a5b49`. The working tree carried uncommitted KiCad/doc
  edits under `hardware/phantasm/` and `docs/` only; nothing on the firmware
  build path was dirty.
- Numbers are stable across the 664 commits since the 2026-07-27 capture at
  `ece0955b`: worst hold 47.79 → 47.05 ms, peak render 56.97 → 56.91 ms,
  spilled 0/3328 either way. The colour-refactor and OpLeg/mesh review batches
  in that range cost nothing measurable here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=IslamicStars`,
`HS_PROFILE_WINDOW=16`; `just profile IslamicStars` = build + flash + capture.
