# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-28, **-O3**)

Point-in-time snapshot (regenerate with the explicit `profile_o3` Reproduce command).
Raw capture: `build/prof/islamicstars_o3.log`; retained raw evidence (supporting artifact removed).
Replaces the prior 2026-09-24 unverified candidate report with a clean committed-source capture.
This standard report measures the current baseline. The unlanded convex-face
AA candidate is measured separately in the matched comparison (supporting artifact removed).

[Shipping sibling](../shipping/profile_islamicstars_teensy_2026-09-28.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile_o3`: global -O3 -ffast-math single-effect reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars, 288×144, single-entry playlist, source `0c02f3912a98677184d3cf43b31a5d468d4b8ba3` |
| Method | 210 s, window 16, `HS_PROFILE_EPOCH_REVS=4000`, `-D HS_PROFILE_TRANS_SPEED=4`; runtime frames 2–3337; scope/ISR windows 17–3328; captured 2026-09-28 16:50 local time |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh IslamicStars profile_o3 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=4000"` |

Image size: `FLASH: code:149296, data:200520, headers:8584` /
`RAM1: variables:315456, code:44872, padding:20664, free:143296` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 2577–2592 root cycles ÷
600 MHz match the measured wall sum within **2.51 ppm**. The untouched
capture passes `tools/parse_profile.py ... validate`, with no epoch reset,
complete per-window frame telemetry, and the expected effect/resolution.

## Frame cadence

README cells: peak 🟢 49.112 (23), spilled 🟢 0/3336 (0.00%).

**Peak live render: 49.112 ms**, frame 2811, marker-owned
shape `truncatedIcosidodecahedron_truncate50d_ambo_dual`. Spilled **0/3336
(0.00%)**; mean render 22.419 ms,
mean wall 62.407 ms. Exact telemetry includes all live
transitions and trailing rows after the final complete scope window.

Setup frame 1 renders in 15.795 ms before publication; it is excluded
from runtime means, peaks, spill counts and denominators. Scope summaries
exclude the whole first window. `is_timeline_step` averages
22.239 ms/frame across the retained complete windows.

One display window is 62.5 ms. The effect renders one 144×72 quadrant,
10,368 pixels, of the 288×144 canvas. `canvas_buffer_wait` is synchronization
idle, not render work. Peak render retains
13.388 ms of margin; all live frames fit the 16 fps budget.

## Phase-by-phase readout

The carousel visits 23 authored recipes. Each shape builds through operator
legs, then holds, ripples, settles and exits. Trans Speed 4 shortens both
holds and build/ripple animation sampling; these are matched TS4 measurements.

### Finished geometry/ripple (frames 2577–2592)

```text
frame                      62.49 ms 37.49 Mcyc 100%
  pov_preserve_half        144.0 us  86.4 kcyc   0% x1 144.0us/c
  is_timeline_step         42.77 ms 25.66 Mcyc  68%
    is_draw_shape          42.72 ms 25.63 Mcyc  68%
      is_mesh_scan         38.71 ms 23.23 Mcyc  62%
        scan_mesh_raster   27.54 ms 16.53 Mcyc  44%
          filter_blend      1.28 ms 769.2 kcyc   2% x18354 41.9cyc/b
        scan_face_setup    10.78 ms  6.47 Mcyc  17% x1082 10.0us/c
      is_face_offsets      527.0 us 316.2 kcyc   1% x1 527.0us/c
      is_mesh_transform     3.47 ms  2.08 Mcyc   6% x1 3474.6us/c
  is_ripple_prepare          9.0 us   5.4 kcyc   0% x1 9.0us/c
  canvas_clear              84.6 us  50.8 kcyc   0% x1 84.6us/c
  canvas_buffer_wait       19.48 ms 11.69 Mcyc  31% x1 19480.2us/c
```

Wall min/avg/max = 60.298/62.491/65.092 ms.
Render averages 43.012 ms. This is the highest render-mean complete
window of this regime; the root includes synchronization idle. Tagged
mixed-parent counters are inclusive shared totals, not exclusive phase costs.

### Build or transition (frames 2561–2576)

```text
frame                      62.59 ms 37.55 Mcyc 100%
  pov_preserve_half        143.5 us  86.1 kcyc   0% x1 143.5us/c
  is_timeline_step         39.19 ms 23.51 Mcyc  63%
    is_build_draw          13.90 ms  8.34 Mcyc  22%
      is_build_scan        13.83 ms  8.30 Mcyc  22% x0 36872.2us/c
      is_mesh_transform     70.4 us  42.2 kcyc   0% x0 187.7us/c
    hk_conway_compile      329.3 us 197.6 kcyc   1% x0 878.1us/c
    hk_conway_sweep        506.6 us 303.9 kcyc   1% x0 1350.9us/c
    is_draw_shape          24.21 ms 14.52 Mcyc  39%
      is_mesh_scan         23.32 ms 13.99 Mcyc  37%
        scan_mesh_raster   25.96 ms 15.58 Mcyc  41%
          filter_blend      1.29 ms 772.0 kcyc   2% x18303 42.2cyc/b
        scan_face_setup    10.79 ms  6.48 Mcyc  17% x1082 10.0us/c
      is_face_offsets      330.3 us 198.2 kcyc   1% x1 528.4us/c
      is_mesh_transform    557.7 us 334.6 kcyc   1% x1 892.4us/c
  is_ripple_prepare          2.9 us   1.7 kcyc   0% x1 2.9us/c
  canvas_clear              84.6 us  50.8 kcyc   0% x1 84.6us/c
  canvas_buffer_wait       23.17 ms 13.90 Mcyc  37% x1 23168.1us/c
```

Wall min/avg/max = 58.809/62.585/65.571 ms.
Render averages 39.418 ms. This is the highest render-mean complete
window of this regime; the root includes synchronization idle. Tagged
mixed-parent counters are inclusive shared totals, not exclusive phase costs.

### Per-preset table

Rows are ranked by the worst clean-hold mesh-scan window. Render and blend
figures come from that same window; a missing blend row is unavailable,
not zero. Windows are clean/owned complete windows. Runtime peak includes
all live frames attributed to the shape, including boundary work.

All 23 recipes appear and wrap at frame 1776. Geometry is the finished
`Built Shape` V/E/F/I, not the spawned seed counts.

Clean holds require no build draw, one owner, no spawn/build-complete marker,
and `scan_mesh_raster` calls equal to 16 × finished face count.

| Shape | V/E/F/I | Windows | Blended px/f | Scan ms | Render ms | Peak ms | fps |
|---|---|---:|---:|---:|---:|---:|---:|
| `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | 3240/4320/1082/8640 | 4/12 | 18,353.8 | 38.713 | 43.012 | 47.666 | 16.00 |
| `truncatedIcosidodecahedron_bevel5_relax_hk77` | 2160/2880/722/5760 | 4/8 | 18,671.5 | 36.379 | 38.436 | 42.339 | 15.89 |
| `truncatedOctahedron_gyro_kis_hk17` | 1620/2160/542/4320 | 4/12 | 18,254.1 | 35.367 | 37.426 | 43.563 | 15.98 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | 1620/2160/542/4320 | 2/8 | 15,997.2 | 35.218 | 37.441 | 39.429 | 16.03 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | 1620/2160/542/4320 | 4/10 | 17,523.9 | 34.220 | 36.282 | 39.769 | 15.95 |
| `truncatedIcosahedron_hk54_ambo_hk72` | 1620/2160/542/4320 | 4/8 | 17,819.0 | 30.339 | 32.443 | 36.290 | 16.06 |
| `truncatedIcosahedron_ambo_relax_truncate33_hk64` | 1620/2160/542/4320 | 4/8 | 16,948.9 | 29.139 | 30.749 | 34.244 | 15.96 |
| `icosahedron_snub_relax_truncate033_hankin62` | 1350/1800/452/3600 | 1/4 | 16,208.1 | 27.090 | 29.060 | 30.761 | 16.04 |
| `dodecahedron_ambo_bevel33_relax_hk66` | 1080/1440/362/2880 | 4/10 | 16,050.8 | 26.193 | 27.823 | 29.760 | 16.01 |
| `truncatedIcosahedron_hk58_chamfer63` | 990/1440/452/2880 | 4/8 | 16,474.1 | 25.925 | 27.449 | 28.550 | 16.01 |
| `dodecahedron_hk72_ambo_dual_hk20` | 540/720/182/1440 | 2/5 | 14,766.6 | 25.201 | 25.939 | 29.844 | 15.91 |
| `rhombicuboctahedron_hk63_ambo_hk63` | 864/1152/290/2304 | 4/10 | 15,156.9 | 24.520 | 25.756 | 27.471 | 15.98 |
| `truncatedIcosidodecahedron_truncate50d_ambo_dual` | 542/1080/540/2160 | 2/10 | 18,949.6 | 24.465 | 25.539 | 49.112 | 16.04 |
| `dodecahedron_hk54_ambo_hk72` | 540/720/182/1440 | 4/10 | 14,657.5 | 22.467 | 23.437 | 26.170 | 16.06 |
| `dodecahedron_hk62_ambo_hk62` | 540/720/182/1440 | 4/9 | 14,133.6 | 21.013 | 21.893 | 23.362 | 15.99 |
| `icosahedron_ambo_truncate033_hankin59` | 540/720/182/1440 | 4/8 | 13,752.1 | 20.117 | 20.855 | 22.792 | 15.99 |
| `octahedron_hk17_ambo_hk73` | 216/288/74/576 | 4/8 | 13,112.7 | 19.192 | 19.701 | 22.173 | 16.02 |
| `snubDodecahedron_truncate5d_ambo_dual` | 452/900/450/1800 | 4/10 | 16,917.9 | 18.777 | 19.548 | 37.054 | 16.01 |
| `octahedron_hk34_ambo_hk72` | 216/288/74/576 | 2/8 | 13,134.5 | 18.711 | 19.271 | 22.186 | 16.04 |
| `truncatedIcosahedron_truncate50d_ambo_dual` | 272/540/270/1080 | 2/5 | 15,886.8 | 18.071 | 18.698 | 38.660 | 16.00 |
| `dodecahedron_bevel2_relax_gyro` | 542/900/360/1800 | 4/12 | 15,842.2 | 17.819 | 18.852 | 37.215 | 16.02 |
| `icosahedron_kis_gyro` | 272/450/180/900 | 2/14 | 14,432.2 | 15.160 | 15.828 | 27.851 | 16.02 |
| `icosidodecahedron_truncate5d_ambo_dual` | 182/360/180/720 | 4/10 | 14,413.9 | 14.098 | 14.548 | 22.360 | 16.14 |

### Per-pixel figures

Window 2577–2592 records 18,353.8 blended pixels/frame,
1.77× quadrant coverage, at 41.91 cycles/blend.
Inclusive `is_mesh_scan` costs 1265.57 cycles per blended pixel.
Shared-counter parenting can hide `filter_blend` on other draw paths; this
window-local figure is not an invented whole-run blend total.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame 0.34/1.57/23.50 us CPU 2.90%
isr_pack          144.0/frame 5.99/7.08/10.63 us CPU 1.63%
isr_dma_submit    144.0/frame 0.58/0.94/10.29 us CPU 0.22%
```

Times are per-call min/weighted-average/max. Pack averages
7.08 µs/call; submit 0.94 µs/call.
The 600-byte image/black-strobe transfer is asynchronous; its 24 MHz SPI
bound including byte framing is 230 µs, not CPU submit time. ISR share
4.74% leaves approximately 59.54 ms foreground CPU per
62.5 ms window. Render already includes ISRs; peak render needs
1.000× speedup to fit the display budget.

## Summary ranking

1. `is_timeline_step` — 22.239 ms/frame, 35.6% of root time (inclusive).
2. `scan_mesh_raster` — 18.055 ms/frame, 28.9% of root time (inclusive).
3. `is_build_scan` — 7.103 ms/frame, 11.4% of root time (inclusive).

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
- Trans Speed 4 changes build/ripple sample counts as well as dwell. The
  comparison holds it fixed and does not claim default-speed peak coverage.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`,
`HS_PROFILE_WINDOW=16`; the explicit `profile_o3` Reproduce command runs locked build/flash/capture.
Use the Setup reproduction command for this complete cycle and its flags.

## Global -O3 vs selective -O3

Baseline mean render: 22.790 ms shipping versus 22.419 ms O3
(1.017×). The configs use different boards. Global O3 minus shipping
image size is +23,904 B FLASH code and +8,240 B ITCM.
