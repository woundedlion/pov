# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HankinSolids`).
Raw capture: `build/prof/hankinsolids_ship.log`; [retained raw evidence](../evidence/face_aa_2026-09-28/before_hankinsolids_ship.txt).
Replaces the prior 2026-08-26 report.
This standard report measures the current baseline. The unlanded convex-face
AA candidate is measured separately in the [matched comparison](../face_aa_2026-09-28.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`: -Os base with selective-O3 mesh transforms, SDF Face setup/distance and scan hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids, 288×144, single-entry playlist, source `0c02f3912a98677184d3cf43b31a5d468d4b8ba3` |
| Method | 260 s, window 16, `HS_PROFILE_EPOCH_REVS=4000`, `-D HS_PROFILE_ORDERED_CYCLE`; runtime frames 2–4138; scope/ISR windows 17–4128; captured 2026-09-28 16:45 local time |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HankinSolids profile 260 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=4000"` |

Image size: `FLASH: code:105880, data:170820, headers:8996` /
`RAM1: variables:315392, code:34952, padding:30584, free:143360` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 2241–2256 root cycles ÷
600 MHz match the measured wall sum within **3.39 ppm**. The untouched
capture passes `tools/parse_profile.py ... validate`, with no epoch reset,
complete per-window frame telemetry, and the expected effect/resolution.

## Frame cadence

**Peak live render: 34.837 ms**, frame 2252, marker-owned
shape `truncatedIcosidodecahedron`. Spilled **0/4137
(0.00%)**; mean render 12.739 ms,
mean wall 62.423 ms. Exact telemetry includes all live
transitions and trailing rows after the final complete scope window.

Setup frame 1 renders in 14.512 ms before publication; it is excluded
from runtime means, peaks, spill counts and denominators. Scope summaries
exclude the whole first window. `hk_timeline_step` averages
12.485 ms/frame across the retained complete windows.

One display window is 62.5 ms. The effect renders one 288×36 quadrant,
10,368 pixels, of the 288×144 canvas. `canvas_buffer_wait` is synchronization
idle, not render work. Peak render retains
27.663 ms of margin; all live frames fit the 16 fps budget.

## Phase-by-phase readout

The ordered tour visits all 18 solids over 30 graph legs, including repeated
visits, settling edges and family bridges. Each node has a 64-frame interlace
sweep followed by a 48-frame morph, with 12 settle frames where required.

### Interlace sweep (frames 2241–2256)

```text
frame                      62.76 ms 37.66 Mcyc 100%
  pov_preserve_half        145.9 us  87.6 kcyc   0% x1 145.9us/c
  hk_timeline_step         24.73 ms 14.84 Mcyc  39%
    hk_draw_mesh           23.97 ms 14.38 Mcyc  38%
      hk_mesh_scan         23.93 ms 14.36 Mcyc  38%
        scan_mesh_raster   21.86 ms 13.11 Mcyc  35%
          filter_blend      1.08 ms 648.2 kcyc   2% x15914 40.7cyc/b
        scan_face_setup     1.94 ms  1.17 Mcyc   3% x182 10.7us/c
      hk_mesh_transform     37.4 us  22.4 kcyc   0% x1 37.4us/c
    hk_update_hankin       709.2 us 425.5 kcyc   1% x1 709.2us/c
  canvas_clear              84.2 us  50.5 kcyc   0% x1 84.2us/c
  canvas_buffer_wait       37.80 ms 22.68 Mcyc  60% x1 37800.3us/c
```

Wall min/avg/max = 54.355/62.762/70.526 ms.
Render averages 24.963 ms. This is the highest render-mean complete
window of this regime; the root includes synchronization idle. Tagged
mixed-parent counters are inclusive shared totals, not exclusive phase costs.

### Morph or boundary (frames 2273–2288)

```text
frame                      61.85 ms 37.11 Mcyc 100%
  pov_preserve_half        143.4 us  86.1 kcyc   0% x1 143.4us/c
  hk_timeline_step         17.68 ms 10.61 Mcyc  29%
    hk_draw_mesh            2.53 ms  1.52 Mcyc   4%
      hk_mesh_scan          2.53 ms  1.52 Mcyc   4% x0 10106.2us/c
      hk_mesh_transform      1.8 us   1.1 kcyc   0% x0 7.3us/c
    hk_conway_compile       12.8 us   7.7 kcyc   0% x0 51.2us/c
    hk_conway_sweep         52.1 us  31.3 kcyc   0% x0 208.3us/c
    hk_draw_mesh           14.29 ms  8.57 Mcyc  23%
      hk_mesh_scan         14.26 ms  8.56 Mcyc  23%
        scan_mesh_raster   15.15 ms  9.09 Mcyc  24%
          filter_blend     837.9 us 502.7 kcyc   1% x12974 38.8cyc/b
        scan_face_setup     1.51 ms 908.7 kcyc   2% x152 10.0us/c
      hk_mesh_transform     28.1 us  16.9 kcyc   0% x1 37.5us/c
    hk_update_hankin       492.3 us 295.4 kcyc   1% x1 656.4us/c
  canvas_clear              84.2 us  50.5 kcyc   0% x1 84.2us/c
  canvas_buffer_wait       43.94 ms 26.36 Mcyc  71% x1 43941.1us/c
```

Wall min/avg/max = 54.508/61.851/63.908 ms.
Render averages 17.910 ms. This is the highest render-mean complete
window of this regime; the root includes synchronization idle. Tagged
mixed-parent counters are inclusive shared totals, not exclusive phase costs.

### Per-preset table

Rows are ranked by the worst clean-hold mesh-scan window. Render and blend
figures come from that same window; a missing blend row is unavailable,
not zero. Windows are clean/owned complete windows. Runtime peak includes
all live frames attributed to the shape, including boundary work.

All 18 nodes and 30 landings close the ordered tour at frame 3522. Initial unmarked frames are tetrahedron.

Clean holds require 16 `hk_update_hankin` calls, one owner and no arrival marker.

| Shape | V/E/F/I | Windows | Blended px/f | Scan ms | Render ms | Peak ms | fps |
|---|---|---:|---:|---:|---:|---:|---:|
| `truncatedIcosidodecahedron` | — | 3/8 | 15,913.6 | 23.934 | 24.963 | 34.837 | 15.93 |
| `snubDodecahedron` | — | 4/8 | 14,331.3 | 20.685 | 21.576 | 26.600 | 15.95 |
| `rhombicosidodecahedron` | — | 3/8 | 13,895.4 | 20.341 | 21.134 | 23.198 | 16.04 |
| `truncatedIcosahedron` | — | 3/7 | 13,816.5 | 19.291 | 19.969 | 22.594 | 15.92 |
| `truncatedDodecahedron` | — | 3/7 | 13,812.4 | 19.242 | 19.925 | 23.192 | 15.93 |
| `truncatedCuboctahedron` | — | 4/9 | 13,277.8 | 18.295 | 18.903 | 23.192 | 16.10 |
| `truncatedCube` | — | 3/7 | 12,527.2 | 18.197 | 18.659 | 22.329 | 15.95 |
| `icosidodecahedron` | — | 10/22 | 12,874.9 | 18.086 | 18.625 | 21.382 | 16.00 |
| `snubCube` | — | 3/8 | 12,857.3 | 17.637 | 18.205 | 23.182 | 16.12 |
| `truncatedOctahedron` | — | 6/14 | 12,542.7 | 17.196 | 17.648 | 19.992 | 16.04 |
| `rhombicuboctahedron` | — | 4/7 | 12,393.1 | 17.089 | 17.604 | 19.850 | 15.98 |
| `dodecahedron` | — | 9/22 | 12,341.6 | 16.230 | 16.665 | 18.732 | 16.00 |
| `truncatedTetrahedron` | — | 6/15 | 11,705.8 | 15.803 | 16.191 | 19.656 | 16.09 |
| `cube` | — | 10/22 | 11,492.9 | 15.296 | 15.659 | 17.418 | 15.98 |
| `icosahedron` | — | 7/14 | 12,155.6 | 15.277 | 15.726 | 18.454 | 15.96 |
| `cuboctahedron` | — | 12/29 | 11,803.6 | 15.096 | 15.508 | 17.758 | 16.00 |
| `octahedron` | — | 12/29 | 11,368.8 | 14.393 | 14.761 | 18.040 | 15.93 |
| `tetrahedron` | — | 10/21 | 11,073.1 | 12.601 | 12.937 | 15.009 | 16.03 |

### Per-pixel figures

Window 2241–2256 records 15,913.6 blended pixels/frame,
1.53× quadrant coverage, at 40.73 cycles/blend.
Inclusive `hk_mesh_scan` costs 902.40 cycles per blended pixel.
Shared-counter parenting can hide `filter_blend` on other draw paths; this
window-local figure is not an invented whole-run blend total.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame 0.52/1.63/17.74 us CPU 3.01%
isr_pack          144.0/frame 6.24/6.95/10.12 us CPU 1.60%
isr_dma_submit    144.0/frame 0.58/0.94/10.23 us CPU 0.22%
```

Times are per-call min/weighted-average/max. Pack averages
6.95 µs/call; submit 0.94 µs/call.
The 600-byte image/black-strobe transfer is asynchronous; its 24 MHz SPI
bound including byte framing is 230 µs, not CPU submit time. ISR share
4.82% leaves approximately 59.48 ms foreground CPU per
62.5 ms window. Render already includes ISRs; peak render needs
1.000× speedup to fit the display budget.

## Summary ranking

1. `hk_timeline_step` — 12.485 ms/frame, 20.0% of root time (inclusive).
2. `scan_mesh_raster` — 11.777 ms/frame, 18.9% of root time (inclusive).
3. `hk_update_hankin` — 0.104 ms/frame, 0.2% of root time (inclusive).

No matched WASM/native timing is used; this is the live device result.

## Caveats

- CYCCNT free-runs, so every scope includes ISR time.
- `filter_blend` is registered under its first caller and can be hidden when
  that parent has zero calls. Its per-pixel instrumentation adds overhead.
- Duplicate-name rows represent individual counters with the same label;
  they are not a combined total for every caller of that label.
- Shipping uses the landed selective-O3 transform/SDF/scan regions; global
  O3 is a single-effect reference. Neither relaxes the shipping memory gates.
- Source provenance attests a clean commit. Each before/after pair uses the
  same board; cross-config comparisons use different boards.
- Every live frame is retained, including geometry builds and transitions.
  One capture per revision/config does not establish a confidence interval.
- Ordered profiling replaces random edge choice with the fixed 30-leg tour;
  it preserves each selected leg's drawing work. The epoch is stretched only.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HankinSolids`,
`HS_PROFILE_WINDOW=16`; `just profile HankinSolids` runs locked build/flash/capture.
Use the Setup reproduction command for this complete cycle and its flags.
