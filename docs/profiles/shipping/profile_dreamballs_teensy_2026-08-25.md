# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-08-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DreamBalls`). Raw capture:
`build/prof/dreamballs_ship.log`. This replaces the 2026-08-09 report, which was
captured against a 5-preset roster and is not comparable row-for-row: the roster
is now 10 presets, and unrelated engine work has since moved blended-pixel counts
by ~22% on the shared presets.

Two captures back this report. The shipping numbers are from `bd38081b4`, which
removed the Möbius warp. An earlier capture at `a4ba1fe3e` (warp still in) is
retained below for the two questions it answered: the cost of the segmented
clip-parameterization fix, and why the warp was worth removing.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: shipping `-Os` plus selective-O3 raster/filter regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, 72 px/arm, clip rows 0–71 |
| Effect | DreamBalls 288×144, single-entry playlist, tip `bd38081b4` |
| Method | `HS_PROFILE`, window 16, 230 s, `HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile 230 16 -D HS_PROFILE_EPOCH_REVS=2000` |

Image size (phantasm, at the warp-era tip): `FLASH: code:490944, data:719872,
headers:8768` / `RAM1: variables:314880, code:191608, padding:5000, free:12800`
/ `RAM2: variables:520064, free:4224`.

Exactness cross-check: frames 2849–2864 match root cycles to wall time within
**1.3 ppm**.

## Frame cadence

A display window is 62.5 ms at 480 RPM. Segment 0 renders a 144×72 quadrant —
the **northern** half canvas, row 0 inclusive — alternating x-half each frame.

Pass aggregate: peak frame render **38.73 ms**, spilled **0/3,648 frames (0%)**.
Every preset holds 16 fps, the worst frame retaining 23.8 ms of margin.
`canvas_buffer_wait` is the intentional wait for the next display flip.

## Phase-by-phase readout

Ten 320-frame presets run in order. The capture contains 11 advance markers,
visits all ten, wraps to preset 1, and shows no epoch reset.

### Peak hold — Disdyakis Triacontahedron (frames 2609–2624)

```
frame                    62.52 ms  37.51 Mcyc  100%
  db_timeline_step       37.11 ms  22.27 Mcyc   59%
    db_draw              37.05 ms  22.23 Mcyc   59%
      db_draw_scene      37.05 ms  22.23 Mcyc   59%
        db_mesh_plot     35.60 ms  21.36 Mcyc   57%  x6.0
          filter_blend    5.63 ms   3.38 Mcyc    9%  x42754
        db_displace     1031.5 us  618.9 kcyc    2%  x6.0
        db_orient        256.9 us  154.2 kcyc    0%  x6.0
  canvas_clear            86.6 us   51.9 kcyc    0%
  canvas_buffer_wait     25.32 ms  15.19 Mcyc   40%
```

Wall min/avg/max = 60.44/62.52/64.57 ms; render avg 37.20, max 38.55.

Rasterization dominates: `db_mesh_plot` is 96% of the timeline step. `db_orient`
is what remains of the old `db_warp_orient` once the Möbius transform is gone —
0.26 ms against the 0.91 ms that scope used to cost.

### Cost is now stationary across the hold

The warp-era capture had this preset swinging from 0.41 ms to 62.58 ms within a
single hold. It no longer does:

| preset 9 | min | p50 | p90 | p99 | max | max/p50 |
|---|--:|--:|--:|--:|--:|--:|
| warp (`a4ba1fe3e`) | 0.41 | 37.55 | 50.76 | 59.34 | **62.58** | 1.67× |
| shipped (`bd38081b4`) | 33.86 | **37.01** | 38.04 | 38.53 | **38.73** | **1.05×** |

The median barely moved (37.55 → 37.01) while the peak fell **38%**. That is the
signature of a redistribution being removed rather than of work being removed:
the warp was not adding average cost, it was concentrating the wireframe into
some frames of some segment bands and vacating others. Per-frame cost across the
hold now varies by 5%.

### Per-preset table

Rows use each preset's worst clean hold; peak is the worst single frame owned by
that preset. Bucket *k* is `PRESETS[k-1]`; bucket 0 is `PRESETS[0]`'s startup
hold before the first marker, which is why 0 and 1 agree. Cycle wrap to preset 1
confirmed, 20 windows per preset (28 for bucket 2).

| # | Solid | copies | render ms | peak ms | fps |
|---|---|--:|--:|--:|--:|
| 9 | Disdyakis Triacontahedron | 6 | 37.29 | 38.73 | 16 |
| 1 | Rhombicuboctahedron | 18 | 35.59 | 36.55 | 16 |
| 0 | Rhombicuboctahedron (startup) | 18 | 35.19 | 36.06 | 16 |
| 10 | Triakis Icosahedron | 6 | 27.30 | 29.49 | 16 |
| 8 | Triakis Icosahedron | 6 | 27.76 | 29.36 | 16 |
| 4 | Icosidodecahedron | 10 | 21.54 | 24.18 | 16 |
| 3 | Truncated Cuboctahedron | 6 | 19.97 | 21.64 | 16 |
| 2 | Rhombicosidodecahedron | 6 | 19.62 | 20.48 | 16 |
| 7 | Triakis Icosahedron | 4.5 | 18.37 | 19.17 | 16 |
| 6 | Truncated Dodecahedron | 4.5 | 12.88 | 14.16 | 16 |
| 5 | Snub Cube | 4.5 | 12.62 | 13.21 | 16 |

Every preset is green, and the spread between a preset's clean hold and its worst
frame is now 1–3 ms rather than the 13.7 ms preset 9 carried with the warp.

### Per-pixel figures

The peak hold blends 42,754 samples/frame against a 10,368-pixel quadrant (4.1×
overdraw) at 79 cyc/blend inside `filter_blend`. Whole-`db_mesh_plot` cost is 500
cyc per blended sample, so the geodesic walk, clipping and shading — not the
blend — set the cadence. Blends per frame are *higher* than the warp-era peak
window's 36,827: without the warp this band sees its even share of the mesh every
frame instead of feast-and-famine.

## Column-ISR / DMA marshaling cost

```
isr_wake        1154/frame  min/avg/max 0.59/1.68/10.96 us  cpu 3.10%
isr_pack         144/frame  min/avg/max 6.23/6.80/9.12 us   cpu 1.56%
isr_dma_submit   144/frame  min/avg/max 0.79/0.95/1.98 us   cpu 0.21%
```

ISR time is included in the foreground scopes. Total ISR CPU share is 4.87%,
leaving ~59.5 ms of the 62.5 ms window for render against a 38.73 ms worst frame.

## Summary ranking

1. `db_mesh_plot` — 35.60 ms in the peak window, 96% of the timeline step; the
   geodesic walk dominates and sets cadence.
2. `filter_blend` — 5.63 ms, already included beneath mesh plotting.
3. `db_displace` — 1.03 ms; `db_orient` adds 0.26 ms.

## Prior capture at `a4ba1fe3e` (warp still in)

### Cost of the segmented clip-parameterization fix

`a4ba1fe3e` stopped `Plot::Mesh::draw_edge` re-parameterizing a clipped edge per
band. A/B at the same tip with only that commit's two files reverted, both
captures window 16 / 230 s / epoch 2000 on the same roster (baseline COM3, fixed
COM4):

| # | Solid | render before → after | Δ | peak before → after |
|---|---|--:|--:|--:|
| 9 | Disdyakis Triacontahedron | 48.81 → 48.90 | +0.2% | 62.37 → 62.58 |
| 0 | Rhombicuboctahedron | 39.19 → 40.31 | +2.9% | 46.70 → 47.92 |
| 1 | Rhombicuboctahedron | 37.85 → 38.95 | +2.9% | 42.61 → 43.62 |
| 4 | Icosidodecahedron | 23.60 → 24.10 | +2.1% | 33.60 → 34.03 |
| 10 | Triakis Icosahedron | 28.04 → 28.41 | +1.3% | 30.49 → 31.16 |
| 8 | Triakis Icosahedron | 27.69 → 28.14 | +1.6% | 29.43 → 29.68 |
| 7 | Triakis Icosahedron | 22.03 → 22.30 | +1.2% | 28.49 → 28.81 |
| 3 | Truncated Cuboctahedron | 22.72 → 22.92 | +0.9% | 27.18 → 27.36 |
| 2 | Rhombicosidodecahedron | 20.42 → 20.70 | +1.4% | 25.05 → 25.27 |
| 6 | Truncated Dodecahedron | 17.43 → 17.49 | +0.3% | 22.47 → 22.49 |
| 5 | Snub Cube | 12.59 → 12.79 | +1.6% | 13.22 → 13.42 |

Mean **+1.4%** for exact segmented equivalence — shading is preserved and only
the adaptive step walk grows. It took spilled frames 0/3,648 → 1/3,648, but
preset 9 peaked at 62.37 ms *before* the change, 0.13 ms under the window: it was
already on the cliff edge and any comparable increment would have tipped it.
Removing the warp has since restored 23.8 ms of margin, so this is no longer
close.

### Why the warp was removed

`db_warp_orient` applied `mobius_gen.transform()` before the global orientation.
Its own cost was 1.94 ms/frame with `db_displace` — 3.1% of the peak frame — but
a Möbius transform on the sphere compresses one region while expanding another,
so it **redistributed** the wireframe across segment bands. Host replay of preset
9's hold under segment 0's clip, plotted samples/frame:

| | mean | p50 | p95 | max | swing |
|---|--:|--:|--:|--:|--:|
| warp on | 8,789 | 10,887 | 15,347 | **17,911** | min ≈ 0 |
| warp off | 11,799 | 11,781 | 12,211 | **12,455** | **1.1×** |

The device capture at `bd38081b4` confirmed the projection: predicted ~42–45 ms,
measured 38.73.

## Caveats

- CYCCNT scopes include ISR time.
- `filter_blend` is parented beneath `db_mesh_plot`; calls approximate blended samples.
- The capture used a 2,000-revolution epoch; this changes dwell capacity, not frame cost.
- Selective-O3: DreamBalls runs through the raster/filter `HS_O3` regions; it has
  no effect-local region.
- The `a4ba1fe3e` A/B baseline was captured on COM3 and its fixed run on COM4.
  Both boards are Teensy 4.0 @ 600 MHz, but no per-board calibration was
  performed, so that ~1.4% mean delta carries an unquantified inter-board
  component. The shipping capture and the warp-era capture it is compared against
  were both on COM4.
- The image-size line predates the warp removal; the removal only shrinks it.
- Palette assignments changed after this capture (`9a4fdd9b5` gave preset 5 its
  own palette). A palette is a LUT lookup, so per-frame cost is unaffected.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, window 16; the
reproduction command above performs the locked build, flash, and capture.
