# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HankinSolids`).
Raw capture: `build/prof/hankinsolids_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_hankinsolids_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses the Hankin mesh builder, face-SDF, and mesh-scan `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh HankinSolids profile 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size:

```text
FLASH: code:111,888, data:165,732, headers:9,100
       free for files:1,744,896
RAM1:  variables:315,520, code:45,736, padding:19,800
       free for local variables:143,232
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2241–2256 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.7 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hk_timeline_step` averages
15.58 ms/frame; its worst window is 32.63 ms/frame
(frames 2241–2256). Peak frame render is
**45.01 ms** (frames 2241–2256);
spilled **0/3328 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 17.49 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2241–2256)

```text
frame                         63.11 ms  37.87 Mcyc  100%
  hk_timeline_step            32.63 ms  19.58 Mcyc   52%
    hk_draw_mesh              31.90 ms  19.14 Mcyc   51%
      hk_mesh_scan            31.87 ms  19.12 Mcyc   51%  x1.0 31871.8us/c
      hk_mesh_transform        31.8 us  19.11 kcyc    0%  x1.0 31.8us/c
    hk_update_hankin          661.8 us 397.10 kcyc    1%  x1.0 661.8us/c
  canvas_clear                 86.2 us  51.75 kcyc    0%  x1.0 86.2us/c
  canvas_buffer_wait          30.39 ms  18.24 Mcyc   48%  x1.0 30391.9us/c
```

Wall min/avg/max = 54.06/63.11/70.77 ms. `hk_timeline_step`
accounts for 51.7% of this measured frame at 32.63 ms/frame.
The complete render is 32.72 ms; `canvas_buffer_wait` contributes
30.39 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 81–96)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  hk_timeline_step             7.85 ms   4.71 Mcyc   13%
    hk_draw_mesh               7.71 ms   4.63 Mcyc   12%
      hk_mesh_scan             7.71 ms   4.62 Mcyc   12%  x1.0 7707.2us/c
      hk_mesh_transform         1.2 us     781 cyc    0%  x1.0 1.2us/c
    hk_conway_compile          15.2 us   9.17 kcyc    0%  x1.0 15.2us/c
    hk_conway_sweep            24.6 us  14.74 kcyc    0%  x1.0 24.6us/c
  canvas_clear                 84.8 us  50.86 kcyc    0%  x1.0 84.8us/c
  canvas_buffer_wait          54.48 ms  32.69 Mcyc   87%  x1.0 54484.8us/c
```

Wall min/avg/max = 61.39/62.42/63.54 ms. `hk_timeline_step`
accounts for 12.6% of this measured frame at 7.85 ms/frame.
The complete render is 7.93 ms; `canvas_buffer_wait` contributes
54.48 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `hk_timeline_step` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `hk_timeline_step` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `truncatedIcosidodecahedron` | truncatedIcosidodecahedron | 8/8 | — | 32.63 | 45.01 | 16.0 |
| 2 | `snubDodecahedron` | snubDodecahedron | 8/8 | 14,277.0 | 29.12 | 36.56 | 16.0 |
| 3 | `rhombicosidodecahedron` | rhombicosidodecahedron | 8/8 | — | 27.50 | 34.64 | 16.0 |
| 4 | `truncatedIcosahedron` | truncatedIcosahedron | 7/7 | — | 25.90 | 31.91 | 16.0 |
| 5 | `truncatedDodecahedron` | truncatedDodecahedron | 7/7 | — | 25.62 | 33.59 | 16.0 |
| 6 | `snubCube` | snubCube | 8/8 | 12,772.9 | 23.67 | 30.24 | 16.0 |
| 7 | `truncatedCuboctahedron` | truncatedCuboctahedron | 8/8 | — | 23.15 | 30.84 | 16.0 |
| 8 | `icosidodecahedron` | icosidodecahedron | 22/22 | — | 22.73 | 27.10 | 16.0 |
| 9 | `rhombicuboctahedron` | rhombicuboctahedron | 7/7 | 12,365.6 | 22.66 | 28.63 | 16.0 |
| 10 | `truncatedCube` | truncatedCube | 7/7 | — | 22.38 | 30.07 | 16.0 |
| 11 | `dodecahedron` | dodecahedron | 22/22 | — | 22.10 | 24.75 | 16.0 |
| 12 | `icosahedron` | icosahedron | 9/9 | 12,318.5 | 21.02 | 25.13 | 16.0 |
| 13 | `truncatedOctahedron` | truncatedOctahedron | 7/7 | — | 20.99 | 24.69 | 16.0 |
| 14 | `octahedron` | octahedron | 14/14 | 11,439.1 | 20.35 | 23.09 | 16.0 |
| 15 | `cube` | cube | 22/22 | — | 19.68 | 23.88 | 16.0 |
| 16 | `cuboctahedron` | cuboctahedron | 21/21 | — | 19.42 | 23.64 | 16.0 |
| 17 | `truncatedTetrahedron` | truncatedTetrahedron | 8/8 | — | 18.77 | 24.93 | 16.0 |
| 18 | `tetrahedron` | tetrahedron | 8/8 | — | 15.39 | 19.66 | 16.0 |
| 19 | `0` | — | 7/7 | — | 15.17 | 19.76 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `hk_timeline_step` uses 1888.3 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.3/f  0.45/1.65/19.58 us  cpu 3.05%
isr_pack          143.9/f  6.24/7.02/9.76 us  cpu 1.62%
isr_dma_submit    143.9/f  0.58/0.94/4.79 us  cpu 0.22%
```

- Pack plus submit costs 1.15 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.88% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `hk_timeline_step` — 25.0% of aggregated root time, 15.58 ms/frame: inclusive measured scope in the live driver.
2. `hk_draw_mesh` — 19.0% of aggregated root time, 11.87 ms/frame: inclusive measured scope in the live driver.
3. `hk_mesh_scan` — 19.0% of aggregated root time, 11.86 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches the Hankin mesh builder, face-SDF, and mesh-scan `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HankinSolids` and
`HS_PROFILE_WINDOW=16`; `just profile HankinSolids` performs the
locked build, flash, capture, and artifact attestation.
