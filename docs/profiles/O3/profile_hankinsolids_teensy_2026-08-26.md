# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_hankinsolids_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh HankinSolids profile_o3 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"`). Raw capture:
`build/prof/hankinsolids_o3.log`, captured 2026-08-26 01:53 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HankinSolids 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 210 s capture, flags `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920` |
| Reproduce | `bash tools/profile_one.sh HankinSolids profile_o3 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"` |

Image size: `FLASH: code:129,640, data:165,800, headers:8,688` / `RAM1: variables:315,552, code:54,488, padding:11,048, free:143,200` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 2241–2256
root counter cycles ÷ 600 MHz match the measured wall sum within **3.4 ppm**
(`tools/parse_profile.py build/prof/hankinsolids_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hk_timeline_step` avg
15.23 ms/f, worst window 31.37 ms/f
(frames 2241–2256),
peak frame render 42.04 ms, spilled 0/3328
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 2241–2256)

```text
frame                          63.11ms  37.86Mcyc 100.0%
  hk_timeline_step             31.37ms  18.82Mcyc  49.7%
    hk_draw_mesh               30.64ms  18.38Mcyc  48.5%
      hk_mesh_scan             30.60ms  18.36Mcyc  48.5% x1.0 30602us/c
      hk_mesh_transform         32.7us   19.6kcyc   0.1% x1.0 33us/c
    hk_update_hankin           671.5us  402.9kcyc   1.1% x1.0 672us/c
  canvas_clear                  86.2us   51.7kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           31.65ms  18.99Mcyc  50.2% x1.0 31648us/c
```

Wall min/avg/max = 55.80/63.10/68.64 ms. The `hk_timeline_step` scope averages 31.37 ms/f in this window, while the exact frame-render peak is 42.04 ms. That is 20.46 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 1713–1728)

```text
frame                          62.39ms  37.44Mcyc 100.0%
  hk_timeline_step              7.75ms   4.65Mcyc  12.4%
    hk_draw_mesh                7.59ms   4.56Mcyc  12.2%
      hk_mesh_scan              7.59ms   4.56Mcyc  12.2% x1.0 7593us/c
      hk_mesh_transform          1.6us     981cyc   0.0% x1.0 2us/c
    hk_conway_compile            9.9us    6.0kcyc   0.0% x1.0 10us/c
    hk_conway_sweep             76.8us   46.1kcyc   0.1% x1.0 77us/c
  canvas_clear                  85.1us   51.1kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           54.56ms  32.73Mcyc  87.4% x1.0 54557us/c
```

Wall min/avg/max = 61.60/62.40/63.08 ms. The `hk_timeline_step` scope averages 7.75 ms/f in this window, while the exact frame-render peak is 8.41 ms. That is 54.09 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `hk_timeline_step` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `truncatedIcosidodecahedron` | truncatedIcosidodecahedron | 8/8 | — | 31.37 | 42.04 | 16.0 |
| `snubDodecahedron` | snubDodecahedron | 8/8 | 14,335.4 | 27.87 | 35.09 | 16.0 |
| `rhombicosidodecahedron` | rhombicosidodecahedron | 8/8 | — | 26.77 | 30.95 | 16.0 |
| `truncatedIcosahedron` | truncatedIcosahedron | 7/7 | — | 24.81 | 31.02 | 16.0 |
| `truncatedDodecahedron` | truncatedDodecahedron | 7/7 | — | 24.19 | 30.84 | 16.0 |
| `icosidodecahedron` | icosidodecahedron | 22/22 | 12,675.0 | 23.26 | 26.30 | 16.0 |
| `snubCube` | snubCube | 8/8 | 12,861.4 | 22.47 | 30.18 | 16.0 |
| `truncatedCuboctahedron` | truncatedCuboctahedron | 8/8 | — | 22.41 | 29.70 | 16.0 |
| `rhombicuboctahedron` | rhombicuboctahedron | 7/7 | 12,569.4 | 22.19 | 26.95 | 16.0 |
| `truncatedCube` | truncatedCube | 7/7 | — | 21.61 | 26.70 | 16.0 |
| `dodecahedron` | dodecahedron | 22/22 | — | 21.20 | 24.17 | 16.0 |
| `truncatedOctahedron` | truncatedOctahedron | 7/7 | — | 20.28 | 23.61 | 16.0 |
| `cuboctahedron` | cuboctahedron | 21/21 | 11,797.2 | 20.22 | 23.31 | 16.0 |
| `icosahedron` | icosahedron | 9/9 | 12,213.8 | 20.07 | 23.68 | 16.0 |
| `octahedron` | octahedron | 14/14 | 11,409.2 | 19.25 | 22.14 | 16.0 |
| `truncatedTetrahedron` | truncatedTetrahedron | 8/8 | — | 18.23 | 23.91 | 16.0 |
| `cube` | cube | 22/22 | — | 18.19 | 22.79 | 16.0 |
| `tetrahedron` | tetrahedron | 8/8 | 11,027.6 | 15.67 | 19.09 | 16.0 |
| `0` | — | 7/7 | — | 14.85 | 18.80 | 16.0 |

201 shape markers cover 18 distinct entries; cycle closure is confirmed by a repeated shape.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1164.4/frame  min/avg/max 0.40/1.56/11.17 us  cpu 2.88%
isr_pack          145.6/frame  min/avg/max 6.28/7.09/9.38 us  cpu 1.63%
isr_dma_submit    145.6/frame  min/avg/max 0.77/0.93/1.01 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.72% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `hk_timeline_step` — 24.4% of measured root time, 15.23 ms/f average.
2. `hk_draw_mesh` — 18.6% of measured root time, 11.58 ms/f average.
3. `hk_mesh_scan` — 18.6% of measured root time, 11.58 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HankinSolids`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh HankinSolids profile_o3 210 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=1920"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 45.01 ms versus
42.04 ms here: global O3 lowers the peak by 2.97 ms (6.6%). O3-image minus shipping-image
size deltas are **FLASH code +17,752 B** and **ITCM code
+8,752 B**. Spill fractions are compared rather than raw counts:
shipping 0/3328 (0%)
versus O3 0/3328 (0%).
