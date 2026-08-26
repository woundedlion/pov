# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `build/prof/islamicstars_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_islamicstars_teensy_2026-07-28.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `IslamicStars::{transform_shape,draw_shape,draw_build_mesh}`, `Scan::Mesh`, and face-SDF regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | IslamicStars 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size:

```text
FLASH: code:128,192, data:193,176, headers:8,360
       free for files:1,701,888
RAM1:  variables:315,488, code:43,240, padding:22,296
       free for local variables:143,264
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2577–2592 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.4 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `is_timeline_step` averages
26.52 ms/frame; its worst window is 56.48 ms/frame
(frames 2577–2592). Peak frame render is
**64.18 ms** (frames 2801–2816);
spilled **2/3328 frames**
(0.1%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak exceeds one display window by 1.68 ms. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2577–2592)

```text
frame                         62.57 ms  37.54 Mcyc  100%
  is_timeline_step            56.48 ms  33.89 Mcyc   90%
    is_draw_shape             56.43 ms  33.86 Mcyc   90%
      is_mesh_scan            52.38 ms  31.43 Mcyc   84%
        scan_mesh_raster      37.88 ms  22.73 Mcyc   61%
          filter_blend         1.43 ms 855.50 kcyc    2%  x18239.4 0.1us/c
        scan_face_setup       14.13 ms   8.48 Mcyc   23%  x1082.0 13.1us/c
      is_face_offsets         521.3 us 312.81 kcyc    1%  x1.0 521.3us/c
      is_mesh_transform        3.53 ms   2.12 Mcyc    6%  x1.0 3530.1us/c
  is_ripple_prepare             5.4 us   3.29 kcyc    0%  x1.0 5.4us/c
  canvas_clear                 89.5 us  53.72 kcyc    0%  x1.0 89.5us/c
  canvas_buffer_wait           6.00 ms   3.60 Mcyc   10%  x1.0 6003.0us/c
```

Wall min/avg/max = 59.99/62.57/66.04 ms. `is_timeline_step`
accounts for 90.3% of this measured frame at 56.48 ms/frame.
The complete render is 56.57 ms; `canvas_buffer_wait` contributes
6.00 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 2689–2704)

```text
  filter_blend                814.3 us 488.59 kcyc    1%  x11595.2 0.1us/c
frame                         62.48 ms  37.49 Mcyc  100%
  is_timeline_step            10.68 ms   6.41 Mcyc   17%
    is_build_draw             10.20 ms   6.12 Mcyc   16%
      is_build_scan           10.20 ms   6.12 Mcyc   16%  x1.0 10198.8us/c
      is_mesh_transform         2.9 us   1.77 kcyc    0%  x1.0 2.9us/c
    hk_conway_compile          24.0 us  14.40 kcyc    0%  x1.0 24.0us/c
    hk_conway_sweep            66.9 us  40.16 kcyc    0%  x1.0 66.9us/c
  is_ripple_prepare             0.2 us     124 cyc    0%  x1.0 0.2us/c
  canvas_clear                 85.2 us  51.14 kcyc    0%  x1.0 85.2us/c
  canvas_buffer_wait          51.71 ms  31.03 Mcyc   83%  x1.0 51710.9us/c
```

Wall min/avg/max = 60.67/62.48/64.23 ms. `is_timeline_step`
accounts for 17.1% of this measured frame at 10.68 ms/frame.
The complete render is 10.76 ms; `canvas_buffer_wait` contributes
51.71 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `is_timeline_step` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `is_timeline_step` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | V=20 E=30 F=12 I=60 | 12/12 | 18,239.4 | 56.48 | 61.87 | 16.0 |
| 2 | `truncatedIcosidodecahedron_bevel5_relax_hk77` | V=120 E=180 F=62 I=360 | 8/8 | 18,755.6 | 47.97 | 56.20 | 16.0 |
| 3 | `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | V=60 E=90 F=32 I=180 | 8/8 | 16,127.0 | 44.56 | 47.01 | 16.0 |
| 4 | `truncatedIcosidodecahedron_truncate50d_ambo_dual` | V=120 E=180 F=62 I=360 | 10/10 | 20,032.7 | 43.05 | 64.18 | 8.0 |
| 5 | `truncatedOctahedron_gyro_kis_hk17` | V=24 E=36 F=14 I=72 | 12/12 | 18,233.1 | 42.09 | 51.88 | 16.0 |
| 6 | `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | V=60 E=90 F=32 I=180 | 10/10 | 17,496.3 | 41.62 | 44.91 | 16.0 |
| 7 | `truncatedIcosahedron_hk54_ambo_hk72` | V=60 E=90 F=32 I=180 | 8/8 | 17,859.9 | 39.20 | 44.58 | 16.0 |
| 8 | `truncatedIcosahedron_ambo_relax_truncate33_hk64` | V=60 E=90 F=32 I=180 | 8/8 | 16,968.9 | 37.14 | 40.94 | 16.0 |
| 9 | `icosahedron_snub_relax_truncate033_hankin62` | V=12 E=30 F=20 I=60 | 4/4 | 16,599.2 | 34.50 | 36.06 | 16.0 |
| 10 | `snubDodecahedron_truncate5d_ambo_dual` | V=60 E=150 F=92 I=300 | 10/10 | 18,445.0 | 34.36 | 45.47 | 16.0 |
| 11 | `truncatedIcosahedron_truncate50d_ambo_dual` | V=60 E=90 F=32 I=180 | 5/5 | 16,628.6 | 34.27 | 47.88 | 16.0 |
| 12 | `truncatedIcosahedron_hk58_chamfer63` | V=60 E=90 F=32 I=180 | 8/8 | 16,570.1 | 33.45 | 34.70 | 16.0 |
| 13 | `dodecahedron_bevel2_relax_gyro` | V=20 E=30 F=12 I=60 | 12/12 | 17,972.9 | 33.25 | 44.81 | 16.0 |
| 14 | `dodecahedron_ambo_bevel33_relax_hk66` | V=20 E=30 F=12 I=60 | 10/10 | 15,972.4 | 32.11 | 34.23 | 16.0 |
| 15 | `rhombicuboctahedron_hk63_ambo_hk63` | V=24 E=48 F=26 I=96 | 10/10 | 15,325.3 | 29.95 | 34.64 | 16.0 |
| 16 | `dodecahedron_hk54_ambo_hk72` | V=20 E=30 F=12 I=60 | 10/10 | 14,633.1 | 27.02 | 31.65 | 16.0 |
| 17 | `icosahedron_kis_gyro` | V=12 E=30 F=20 I=60 | 14/14 | 16,436.3 | 26.99 | 35.42 | 16.0 |
| 18 | `dodecahedron_hk72_ambo_dual_hk20` | V=20 E=30 F=12 I=60 | 5/5 | 14,848.2 | 26.93 | 33.72 | 16.0 |
| 19 | `dodecahedron_hk62_ambo_hk62` | V=20 E=30 F=12 I=60 | 10/10 | 14,171.7 | 24.56 | 30.94 | 16.0 |
| 20 | `icosahedron_ambo_truncate033_hankin59` | V=12 E=30 F=20 I=60 | 8/8 | 13,969.6 | 23.88 | 27.57 | 16.0 |
| 21 | `icosidodecahedron_truncate5d_ambo_dual` | V=30 E=60 F=32 I=120 | 10/10 | 14,701.7 | 22.37 | 26.74 | 16.0 |
| 22 | `octahedron_hk17_ambo_hk73` | V=6 E=12 F=8 I=24 | 8/8 | 13,223.2 | 21.48 | 23.65 | 16.0 |
| 23 | `octahedron_hk34_ambo_hk72` | V=6 E=12 F=8 I=24 | 8/8 | 13,040.0 | 20.04 | 21.66 | 16.0 |

### Per-pixel figures

The peak dominant-scope window blends 18,239.4 pixels/frame, 1.76× the 10,368-pixel quadrant. `filter_blend` costs 46.9 cycles/blend; `is_timeline_step` uses 1857.8 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.9/f  0.46/1.73/17.69 us  cpu 3.19%
isr_pack          144.0/f  6.24/7.20/12.78 us  cpu 1.66%
isr_dma_submit    144.0/f  0.59/0.94/10.22 us  cpu 0.22%
```

- Pack plus submit costs 1.17 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.07% of one 62.5 ms window. The peak
  render requires 1.03× the one-window budget, so the cadence target
  needs a 1.03× speedup.

## Summary ranking

1. `is_timeline_step` — 42.5% of aggregated root time, 26.52 ms/frame: inclusive measured scope in the live driver.
2. `is_draw_shape` — 27.9% of aggregated root time, 17.40 ms/frame: inclusive measured scope in the live driver.
3. `is_mesh_scan` — 27.2% of aggregated root time, 16.97 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches `IslamicStars::{transform_shape,draw_shape,draw_build_mesh}`, `Scan::Mesh`, and face-SDF regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=1920`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=IslamicStars` and
`HS_PROFILE_WINDOW=16`; `just profile IslamicStars` performs the
locked build, flash, capture, and artifact attestation.
