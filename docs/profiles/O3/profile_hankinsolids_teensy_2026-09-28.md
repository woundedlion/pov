# HankinSolids on-device profile — Teensy 4.0, segmented mode (2026-09-28, **-O3**)

Point-in-time snapshot (regenerate with the explicit `profile_o3` Reproduce command).
Raw capture: `build/prof/hankinsolids_o3.log`; retained raw evidence (supporting artifact removed).
Replaces the prior 2026-08-26 report.
This report measures the September 28 source baseline. The unlanded convex-face
AA candidate is measured separately in the matched comparison (supporting artifact removed).

[Shipping sibling](../shipping/profile_hankinsolids_teensy_2026-10-07.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile_o3`: global -O3 -ffast-math single-effect reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HankinSolids, 288×144, single-entry playlist, source `0c02f3912a98677184d3cf43b31a5d468d4b8ba3` |
| Method | 260 s, window 16, `HS_PROFILE_EPOCH_REVS=4000`, `-D HS_PROFILE_ORDERED_CYCLE`; runtime frames 2–4138; scope/ISR windows 17–4128; captured 2026-09-28 16:45 local time |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HankinSolids profile_o3 260 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=4000"` |

Image size: `FLASH: code:122968, data:170664, headers:8448` /
`RAM1: variables:315392, code:37096, padding:28440, free:143360` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 2241–2256 root cycles ÷
600 MHz match the measured wall sum within **3.25 ppm**. The untouched
capture passes `tools/parse_profile.py ... validate`, with no epoch reset,
complete per-window frame telemetry, and the expected effect/resolution.

## Frame cadence

README cells: peak 🟢 32.865 (18), spilled 🟢 0/4137 (0.00%).

**Peak live render: 32.865 ms**, frame 2252, marker-owned
shape `truncatedIcosidodecahedron`. Spilled **0/4137
(0.00%)**; mean render 12.823 ms,
mean wall 62.429 ms. Exact telemetry includes all live
transitions and trailing rows after the final complete scope window.

Setup frame 1 renders in 14.821 ms before publication; it is excluded
from runtime means, peaks, spill counts and denominators. Scope summaries
exclude the whole first window. `hk_timeline_step` averages
12.570 ms/frame across the retained complete windows.

One display window is 62.5 ms. The effect renders one 144×72 quadrant,
10,368 pixels, of the 288×144 canvas. `canvas_buffer_wait` is synchronization
idle, not render work. Peak render retains
29.635 ms of margin; all live frames fit the 16 fps budget.

## Phase-by-phase readout

The ordered tour visits all 18 solids over 30 graph legs, including repeated
visits, settling edges and family bridges. Each node has a 64-frame interlace
sweep followed by a 48-frame morph, with 12 settle frames where required.

### Interlace sweep (frames 2241–2256)

```text
frame                      62.78 ms 37.67 Mcyc 100%
  pov_preserve_half        145.5 us  87.3 kcyc   0% x1 145.5us/c
  hk_timeline_step         24.41 ms 14.65 Mcyc  39%
    hk_draw_mesh           23.65 ms 14.19 Mcyc  38%
      hk_mesh_scan         23.59 ms 14.15 Mcyc  38%
        scan_mesh_raster   21.52 ms 12.91 Mcyc  34%
          filter_blend      1.07 ms 643.9 kcyc   2% x15836 40.7cyc/b
        scan_face_setup     1.95 ms  1.17 Mcyc   3% x182 10.7us/c
      hk_mesh_transform     57.2 us  34.3 kcyc   0% x1 57.2us/c
    hk_update_hankin       694.5 us 416.7 kcyc   1% x1 694.5us/c
  canvas_clear              84.3 us  50.6 kcyc   0% x1 84.3us/c
  canvas_buffer_wait       38.14 ms 22.88 Mcyc  61% x1 38139.4us/c
```

Wall min/avg/max = 55.767/62.778/69.095 ms.
Render averages 24.640 ms. This is the highest render-mean complete
window of this regime; the root includes synchronization idle. Tagged
mixed-parent counters are inclusive shared totals, not exclusive phase costs.

### Morph or boundary (frames 1217–1232)

```text
frame                      62.05 ms 37.23 Mcyc 100%
  pov_preserve_half        142.1 us  85.3 kcyc   0% x1 142.1us/c
  hk_timeline_step         17.80 ms 10.68 Mcyc  29%
    hk_draw_mesh           582.3 us 349.4 kcyc   1%
      hk_mesh_scan         578.8 us 347.3 kcyc   1% x0 9260.1us/c
      hk_mesh_transform      1.4 us   0.8 kcyc   0% x0 21.7us/c
    hk_conway_compile        2.8 us   1.7 kcyc   0% x0 44.7us/c
    hk_conway_sweep          7.3 us   4.4 kcyc   0% x0 116.4us/c
    hk_draw_mesh           16.79 ms 10.08 Mcyc  27%
      hk_mesh_scan         16.76 ms 10.06 Mcyc  27%
        scan_mesh_raster   16.61 ms  9.97 Mcyc  27%
          filter_blend     833.5 us 500.1 kcyc   1% x12865 38.9cyc/b
        scan_face_setup    644.6 us 386.7 kcyc   1% x60 10.7us/c
      hk_mesh_transform     25.6 us  15.3 kcyc   0% x1 27.3us/c
    hk_update_hankin       225.5 us 135.3 kcyc   0% x1 240.5us/c
  canvas_clear              84.3 us  50.6 kcyc   0% x1 84.3us/c
  canvas_buffer_wait       44.02 ms 26.41 Mcyc  71% x1 44021.1us/c
```

Wall min/avg/max = 54.017/62.047/64.657 ms.
Render averages 18.026 ms. This is the highest render-mean complete
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
| `truncatedIcosidodecahedron` | — | 3/8 | 15,835.9 | 23.589 | 24.640 | 32.865 | 15.93 |
| `snubDodecahedron` | — | 4/8 | 14,250.6 | 20.602 | 21.535 | 26.173 | 15.93 |
| `rhombicosidodecahedron` | — | 3/8 | 13,836.2 | 19.879 | 20.696 | 22.968 | 16.03 |
| `truncatedIcosahedron` | — | 3/7 | 13,841.8 | 19.372 | 20.072 | 22.533 | 15.92 |
| `truncatedDodecahedron` | — | 3/7 | 13,746.8 | 19.156 | 19.855 | 23.749 | 15.94 |
| `truncatedCuboctahedron` | — | 4/9 | 13,302.0 | 18.293 | 18.944 | 23.175 | 16.10 |
| `truncatedCube` | — | 3/7 | 12,466.4 | 18.034 | 18.516 | 21.712 | 15.94 |
| `icosidodecahedron` | — | 10/22 | 12,798.3 | 17.596 | 18.180 | 20.399 | 16.01 |
| `snubCube` | — | 3/8 | 12,823.2 | 17.517 | 18.115 | 23.052 | 16.12 |
| `truncatedOctahedron` | — | 6/14 | 12,560.2 | 17.313 | 17.792 | 19.554 | 16.04 |
| `rhombicuboctahedron` | — | 4/7 | 12,327.9 | 16.808 | 17.337 | 19.616 | 15.98 |
| `truncatedTetrahedron` | — | 6/15 | 11,755.6 | 16.135 | 16.559 | 19.718 | 16.06 |
| `dodecahedron` | — | 9/22 | 12,296.7 | 15.879 | 16.346 | 18.292 | 15.99 |
| `icosahedron` | — | 7/14 | 12,167.1 | 15.699 | 16.162 | 18.757 | 15.94 |
| `cuboctahedron` | — | 12/29 | 11,765.5 | 15.542 | 15.989 | 17.869 | 16.02 |
| `cube` | — | 10/22 | 11,431.9 | 15.138 | 15.532 | 17.814 | 15.99 |
| `octahedron` | — | 12/29 | 11,432.8 | 14.829 | 15.230 | 17.952 | 16.06 |
| `tetrahedron` | — | 10/21 | 11,059.6 | 12.353 | 12.721 | 14.557 | 16.02 |

### Per-pixel figures

Window 2241–2256 records 15,835.9 blended pixels/frame,
1.53× quadrant coverage, at 40.66 cycles/blend.
Inclusive `hk_mesh_scan` costs 893.75 cycles per blended pixel.
Shared-counter parenting can hide `filter_blend` on other draw paths; this
window-local figure is not an invented whole-run blend total.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame 0.37/1.49/30.62 us CPU 2.74%
isr_pack          144.0/frame 5.99/6.86/9.85 us CPU 1.58%
isr_dma_submit    144.0/frame 0.59/0.93/11.26 us CPU 0.21%
```

Times are per-call min/weighted-average/max. Pack averages
6.86 µs/call; submit 0.93 µs/call.
The 600-byte image/black-strobe transfer is asynchronous; its 24 MHz SPI
bound including byte framing is 230 µs, not CPU submit time. Inclusive
`isr_wake` share 2.74% (pack and submit run nested inside it)
leaves approximately 60.79 ms foreground CPU per
62.5 ms window before unmeasured interrupts. Render already includes ISRs; peak render needs
1.000× speedup to fit the display budget.

## Summary ranking

1. `hk_timeline_step` — 12.570 ms/frame, 20.1% of root time (inclusive).
2. `scan_mesh_raster` — 11.853 ms/frame, 19.0% of root time (inclusive).
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
`HS_PROFILE_WINDOW=16`; the explicit `profile_o3` Reproduce command runs locked build/flash/capture.
Use the Setup reproduction command for this complete cycle and its flags.

## Global -O3 vs selective -O3

Baseline mean render: 12.739 ms shipping versus 12.823 ms O3
(0.993×). The configs use different boards. Global O3 minus shipping
image size is +17,088 B FLASH code and +2,144 B ITCM.

This historical global-O3 comparison and its image deltas retain the September 28 source pair at `0c02f3912a98677184d3cf43b31a5d468d4b8ba3`. The shipping sibling link points to the later October 7 sector-distance snapshot.
