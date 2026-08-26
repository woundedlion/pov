# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `build/prof/islamicstars_ship.log`, sourced from the isolated
regression-reclaim tree. Replaces the earlier same-day sweep capture after
the cadence-reclamation changes in `0df961b81`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `IslamicStars::{transform_shape,draw_shape,draw_build_mesh}`, `Scan::Mesh`, and face-SDF regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `0df961b818ae08c5f58639edb93d812164a9356a` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size:

```text
FLASH: code:128,192, data:193,176, headers:8,360
       free for files:1,701,888
RAM1:  variables:315,488, code:43,496, padding:22,040
       free for local variables:143,264
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2577–2592 root counter cycles ÷ 600 MHz
match the measured wall sum within **4.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, a repeated shape closing the cycle, and no epoch
reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`is_timeline_step` averages 22.21 ms/frame; its worst window is
43.22 ms/frame (frames 2577–2592). Peak frame render is **50.65 ms**
(frames 2801–2816); spilled **0/3328 frames (0.0%)**.

A display window is 62.5 ms. The segmented master renders one quadrant,
approximately 10,368 pixels. Every captured shape and transition holds 16 fps;
the peak retains 11.85 ms of the display window. The `canvas_buffer_wait`
scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the compressed carousel covers 23 distinct shapes and repeats
its first shape, closing the cycle. The sections below bound its held-draw and
build regimes; per-shape ownership includes the following transition.

### Highest-cost clean hold (window frames 2577–2592)

```text
frame                        62.56ms 37.54Mcyc 100%
  is_timeline_step           43.22ms 25.93Mcyc  69%
    is_draw_shape            43.17ms 25.90Mcyc  69%
      is_mesh_scan           39.12ms 23.47Mcyc  63%
        scan_mesh_raster     25.51ms 15.30Mcyc  41%
          filter_blend        1.31ms 788.94kcyc  2% x18214 43.3cyc/c
        scan_face_setup      13.22ms  7.93Mcyc  21% x1082 12.2us/c
      is_face_offsets        522.3us 313.36kcyc   1% x1.0 522.3us/c
      is_mesh_transform       3.53ms  2.12Mcyc   6% x1.0 3529.9us/c
  is_ripple_prepare            5.8us  3.46kcyc   0% x1.0 5.8us/c
  canvas_clear                87.4us 52.44kcyc   0% x1.0 87.4us/c
  canvas_buffer_wait         19.24ms 11.54Mcyc  31% x1.0 19240.3us/c
```

Wall min/avg/max = 60.55/62.56/64.52 ms. `is_timeline_step` accounts for
69.1% of this measured frame at 43.22 ms/frame. Complete render averages
43.32 ms and peaks at 48.97 ms for this shape; `canvas_buffer_wait` contributes
19.24 ms of flip-alignment idle. The pass-wide 50.65 ms peak belongs to the
transition owned by `truncatedIcosidodecahedron_truncate50d_ambo_dual`.

### Lowest-cost build regime (window frames 2689–2704)

```text
  filter_blend               738.6us 443.21kcyc   1% x11586 38.3cyc/c
frame                        62.49ms 37.49Mcyc 100%
  is_timeline_step            9.31ms  5.58Mcyc  15%
    is_build_draw             8.83ms  5.30Mcyc  14%
      is_build_scan           8.83ms  5.30Mcyc  14% x1.0 8831.4us/c
      is_mesh_transform        2.9us  1.73kcyc   0% x1.0 2.9us/c
    hk_conway_compile         21.8us 13.07kcyc   0% x1.0 21.8us/c
    hk_conway_sweep           65.4us 39.27kcyc   0% x1.0 65.4us/c
  is_ripple_prepare            0.2us    122cyc   0% x1.0 0.2us/c
  canvas_clear                85.4us 51.28kcyc   0% x1.0 85.4us/c
  canvas_buffer_wait         53.09ms 31.86Mcyc  85% x1.0 53092.6us/c
```

Wall min/avg/max = 60.23/62.49/64.52 ms. This build leg spends 9.31 ms/frame
in the timeline and 8.83 ms in `is_build_scan`; complete render averages
9.39 ms. The top-level `filter_blend` row is the documented first-entry
parenting artifact, not work outside the frame.

### Per-preset table

The parser observed a closed cycle: the marker sequence repeated its first
shape. Rows rank modal-call-count clean holds by `is_timeline_step` cost;
geometry comes from each `Spawning Shape` marker. Cadence peaks use per-frame
ownership, including the following transition. All 208 owned windows meet the
modal clean-hold criterion.

| rank | shape | V/E/F/I | windows | blended px/f | `is_timeline_step` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | 20/30/12/60 | 12/12 | 18,214.2 | 43.22 | 48.97 | 16.0 |
| 2 | `truncatedIcosidodecahedron_bevel5_relax_hk77` | 120/180/62/360 | 8/8 | 18,755.6 | 39.61 | 43.42 | 16.0 |
| 3 | `truncatedIcosidodecahedron_truncate50d_ambo_dual` | 120/180/62/360 | 10/10 | 20,032.7 | 37.73 | 50.65 | 16.0 |
| 4 | `truncatedOctahedron_gyro_kis_hk17` | 24/36/14/72 | 12/12 | 18,233.1 | 36.34 | 44.33 | 16.0 |
| 5 | `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | 60/90/32/180 | 8/8 | 16,127.0 | 35.96 | 39.02 | 16.0 |
| 6 | `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | 60/90/32/180 | 10/10 | 17,496.3 | 35.03 | 38.15 | 16.0 |
| 7 | `truncatedIcosahedron_hk54_ambo_hk72` | 60/90/32/180 | 8/8 | 17,859.9 | 32.40 | 37.45 | 16.0 |
| 8 | `snubDodecahedron_truncate5d_ambo_dual` | 60/150/92/300 | 10/10 | 18,444.1 | 30.58 | 41.33 | 16.0 |
| 9 | `truncatedIcosahedron_ambo_relax_truncate33_hk64` | 60/90/32/180 | 8/8 | 16,978.2 | 30.48 | 33.59 | 16.0 |
| 10 | `dodecahedron_bevel2_relax_gyro` | 20/30/12/60 | 12/12 | 17,972.9 | 29.32 | 40.94 | 16.0 |
| 11 | `icosahedron_snub_relax_truncate033_hankin62` | 12/30/20/60 | 4/4 | 16,530.4 | 28.36 | 29.86 | 16.0 |
| 12 | `truncatedIcosahedron_hk58_chamfer63` | 60/90/32/180 | 8/8 | 16,558.1 | 27.94 | 29.23 | 16.0 |
| 13 | `truncatedIcosahedron_truncate50d_ambo_dual` | 60/90/32/180 | 5/5 | 16,575.6 | 27.84 | 39.75 | 16.0 |
| 14 | `dodecahedron_ambo_bevel33_relax_hk66` | 20/30/12/60 | 10/10 | 15,972.4 | 26.56 | 28.75 | 16.0 |
| 15 | `rhombicuboctahedron_hk63_ambo_hk63` | 24/48/26/96 | 10/10 | 15,325.3 | 25.17 | 28.71 | 16.0 |
| 16 | `icosahedron_kis_gyro` | 12/30/20/60 | 14/14 | 16,436.3 | 23.59 | 28.92 | 16.0 |
| 17 | `dodecahedron_hk72_ambo_dual_hk20` | 20/30/12/60 | 5/5 | 14,697.9 | 23.00 | 28.20 | 16.0 |
| 18 | `dodecahedron_hk54_ambo_hk72` | 20/30/12/60 | 10/10 | 14,633.1 | 22.94 | 26.94 | 16.0 |
| 19 | `dodecahedron_hk62_ambo_hk62` | 20/30/12/60 | 10/10 | 14,171.7 | 20.72 | 22.44 | 16.0 |
| 20 | `icosahedron_ambo_truncate033_hankin59` | 12/30/20/60 | 8/8 | 13,951.7 | 19.91 | 23.73 | 16.0 |
| 21 | `icosidodecahedron_truncate5d_ambo_dual` | 30/60/32/120 | 10/10 | 14,701.7 | 19.78 | 23.31 | 16.0 |
| 22 | `octahedron_hk17_ambo_hk73` | 6/12/8/24 | 8/8 | 13,223.2 | 18.34 | 20.05 | 16.0 |
| 23 | `octahedron_hk34_ambo_hk72` | 6/12/8/24 | 8/8 | 13,040.0 | 17.08 | 18.35 | 16.0 |

### Per-pixel figures

The highest-cost clean hold blends 18,214.2 pixels/frame, 1.76× the 10,368
pixel quadrant. `filter_blend` costs 43.3 cycles/blend, while
`scan_mesh_raster` uses 840.2 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.2/f  0.45/1.72/23.00us  cpu 3.18%
isr_pack          143.9/f  6.23/7.17/11.89us  cpu 1.65%
isr_dma_submit    143.9/f  0.59/0.94/7.97us   cpu 0.22%
```

- Pack plus submit costs 1.17 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling
  only.
- Aggregate ISR CPU share is 5.05% of one 62.5 ms window. The peak render
  requires 0.81× the one-window budget, so the cadence target needs no
  speedup.

## Summary ranking

1. `is_timeline_step` — 35.6% of aggregate root time, 22.21 ms/frame.
2. `is_draw_shape` — 23.1% of aggregate root time, 14.40 ms/frame.
3. `is_mesh_scan` — 22.4% of aggregate root time, 13.97 ms/frame.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches
  `IslamicStars::{transform_shape,draw_shape,draw_build_mesh}`, `Scan::Mesh`,
  and face-SDF `HS_O3` regions; the rest retains the `-Os` base policy.
- `HS_PROFILE_TRANS_SPEED=4` compresses dwell/transition time and
  `HS_PROFILE_EPOCH_REVS=1920` keeps the full cycle in one epoch. Neither knob
  changes per-frame shape cost.
- Provenance attests a clean source tree at
  `0df961b818ae08c5f58639edb93d812164a9356a`; no uncommitted source state was
  profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=IslamicStars` and
`HS_PROFILE_WINDOW=16`; `just profile IslamicStars` performs the locked build,
flash, capture, and artifact attestation.
