# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `build/prof/islamicstars_ship.log`; [retained raw evidence](../evidence/face_aa_2026-09-28/before_islamicstars_ship.txt).
Replaces the prior 2026-09-24 unverified candidate report with a clean committed-source capture.
This standard report measures the current baseline. The unlanded convex-face
AA candidate is measured separately in the [matched comparison](../face_aa_2026-09-28.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`: -Os base with selective-O3 mesh transforms, SDF Face setup/distance and scan hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars, 288×144, single-entry playlist, source `0c02f3912a98677184d3cf43b31a5d468d4b8ba3` |
| Method | 210 s, window 16, `HS_PROFILE_EPOCH_REVS=4000`, `-D HS_PROFILE_TRANS_SPEED=4`; runtime frames 2–3337; scope/ISR windows 17–3328; captured 2026-09-28 16:50 local time |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh IslamicStars profile 210 16 "-D HS_PROFILE_TRANS_SPEED=4 -D HS_PROFILE_EPOCH_REVS=4000"` |

Image size: `FLASH: code:125392, data:201084, headers:8372` /
`RAM1: variables:315424, code:36632, padding:28904, free:143328` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 2577–2592 root cycles ÷
600 MHz match the measured wall sum within **2.36 ppm**. The untouched
capture passes `tools/parse_profile.py ... validate`, with no epoch reset,
complete per-window frame telemetry, and the expected effect/resolution.

## Frame cadence

**Peak live render: 50.500 ms**, frame 2811, marker-owned
shape `truncatedIcosidodecahedron_truncate50d_ambo_dual`. Spilled **0/3336
(0.00%)**; mean render 22.790 ms,
mean wall 62.423 ms. Exact telemetry includes all live
transitions and trailing rows after the final complete scope window.

Setup frame 1 renders in 15.628 ms before publication; it is excluded
from runtime means, peaks, spill counts and denominators. Scope summaries
exclude the whole first window. `is_timeline_step` averages
22.615 ms/frame across the retained complete windows.

One display window is 62.5 ms. The effect renders one 288×36 quadrant,
10,368 pixels, of the 288×144 canvas. `canvas_buffer_wait` is synchronization
idle, not render work. Peak render retains
12.000 ms of margin; all live frames fit the 16 fps budget.

## Phase-by-phase readout

The carousel visits 23 authored recipes. Each shape builds through operator
legs, then holds, ripples, settles and exits. Trans Speed 4 shortens both
holds and build/ripple animation sampling; these are matched TS4 measurements.

### Finished geometry/ripple (frames 2577–2592)

```text
frame                      62.36 ms 37.42 Mcyc 100%
  pov_preserve_half        142.1 us  85.3 kcyc   0% x1 142.1us/c
  is_timeline_step         44.69 ms 26.81 Mcyc  72%
    is_draw_shape          44.62 ms 26.77 Mcyc  72%
      is_mesh_scan         40.58 ms 24.35 Mcyc  65%
        scan_mesh_raster   28.84 ms 17.30 Mcyc  46%
          filter_blend      1.29 ms 772.7 kcyc   2% x18491 41.8cyc/b
        scan_face_setup    11.33 ms  6.80 Mcyc  18% x1082 10.5us/c
      is_face_offsets      517.3 us 310.4 kcyc   1% x1 517.3us/c
      is_mesh_transform     3.52 ms  2.11 Mcyc   6% x1 3522.4us/c
  is_ripple_prepare          7.1 us   4.2 kcyc   0% x1 7.1us/c
  canvas_clear              84.7 us  50.8 kcyc   0% x1 84.7us/c
  canvas_buffer_wait       17.44 ms 10.46 Mcyc  28% x1 17437.7us/c
```

Wall min/avg/max = 58.903/62.358/67.667 ms.
Render averages 44.922 ms. This is the highest render-mean complete
window of this regime; the root includes synchronization idle. Tagged
mixed-parent counters are inclusive shared totals, not exclusive phase costs.

### Build or transition (frames 2561–2576)

```text
frame                      62.63 ms 37.58 Mcyc 100%
  pov_preserve_half        141.2 us  84.7 kcyc   0% x1 141.2us/c
  is_timeline_step         40.20 ms 24.12 Mcyc  64%
    is_build_draw          14.15 ms  8.49 Mcyc  23%
      is_build_scan        14.08 ms  8.45 Mcyc  22% x0 37539.2us/c
      is_mesh_transform     72.6 us  43.6 kcyc   0% x0 193.7us/c
    hk_conway_compile      357.3 us 214.4 kcyc   1% x0 952.7us/c
    hk_conway_sweep        571.3 us 342.8 kcyc   1% x0 1523.4us/c
    is_draw_shape          24.85 ms 14.91 Mcyc  40%
      is_mesh_scan         23.97 ms 14.38 Mcyc  38%
        scan_mesh_raster   26.44 ms 15.87 Mcyc  42%
          filter_blend      1.28 ms 768.9 kcyc   2% x18265 42.1cyc/b
        scan_face_setup    11.20 ms  6.72 Mcyc  18% x1082 10.3us/c
      is_face_offsets      325.5 us 195.3 kcyc   1% x1 520.8us/c
      is_mesh_transform    559.5 us 335.7 kcyc   1% x1 895.3us/c
  is_ripple_prepare          2.8 us   1.7 kcyc   0% x1 2.8us/c
  canvas_clear              84.7 us  50.8 kcyc   0% x1 84.7us/c
  canvas_buffer_wait       22.19 ms 13.32 Mcyc  35% x1 22193.6us/c
```

Wall min/avg/max = 59.267/62.627/66.059 ms.
Render averages 40.434 ms. This is the highest render-mean complete
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
| `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` | 3240/4320/1082/8640 | 4/12 | 18,490.8 | 40.585 | 44.922 | 49.443 | 16.04 |
| `truncatedIcosidodecahedron_bevel5_relax_hk77` | 2160/2880/722/5760 | 4/8 | 18,671.3 | 37.063 | 39.141 | 43.097 | 15.88 |
| `truncatedOctahedron_gyro_kis_hk17` | 1620/2160/542/4320 | 4/12 | 18,253.2 | 35.887 | 37.972 | 44.549 | 15.98 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin59` | 1620/2160/542/4320 | 2/8 | 15,946.7 | 35.675 | 37.929 | 41.868 | 16.11 |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin73` | 1620/2160/542/4320 | 4/10 | 17,519.3 | 34.666 | 36.749 | 40.246 | 15.95 |
| `truncatedIcosahedron_hk54_ambo_hk72` | 1620/2160/542/4320 | 4/8 | 17,771.8 | 31.260 | 33.395 | 38.062 | 16.07 |
| `truncatedIcosahedron_ambo_relax_truncate33_hk64` | 1620/2160/542/4320 | 4/8 | 16,659.3 | 29.214 | 30.847 | 33.075 | 15.95 |
| `icosahedron_snub_relax_truncate033_hankin62` | 1350/1800/452/3600 | 1/4 | 16,488.2 | 27.786 | 29.781 | 31.150 | 16.04 |
| `truncatedIcosahedron_hk58_chamfer63` | 990/1440/452/2880 | 4/8 | 16,489.2 | 26.545 | 28.092 | 29.832 | 16.02 |
| `dodecahedron_ambo_bevel33_relax_hk66` | 1080/1440/362/2880 | 4/10 | 16,026.6 | 26.068 | 27.718 | 30.080 | 15.99 |
| `truncatedIcosidodecahedron_truncate50d_ambo_dual` | 542/1080/540/2160 | 2/10 | 18,883.2 | 25.423 | 26.520 | 50.500 | 16.06 |
| `rhombicuboctahedron_hk63_ambo_hk63` | 864/1152/290/2304 | 4/10 | 15,132.6 | 24.213 | 25.468 | 27.568 | 15.97 |
| `dodecahedron_hk72_ambo_dual_hk20` | 540/720/182/1440 | 2/5 | 14,711.4 | 23.830 | 24.680 | 29.701 | 16.02 |
| `dodecahedron_hk54_ambo_hk72` | 540/720/182/1440 | 4/10 | 14,646.4 | 21.829 | 22.803 | 23.378 | 16.00 |
| `dodecahedron_hk62_ambo_hk62` | 540/720/182/1440 | 4/9 | 14,125.8 | 21.143 | 22.039 | 23.193 | 15.99 |
| `icosahedron_ambo_truncate033_hankin59` | 540/720/182/1440 | 4/8 | 13,943.4 | 20.705 | 21.462 | 23.258 | 15.97 |
| `snubDodecahedron_truncate5d_ambo_dual` | 452/900/450/1800 | 4/10 | 16,833.3 | 19.108 | 19.979 | 40.601 | 16.02 |
| `octahedron_hk17_ambo_hk73` | 216/288/74/576 | 4/8 | 13,070.6 | 19.036 | 19.537 | 21.643 | 16.01 |
| `octahedron_hk34_ambo_hk72` | 216/288/74/576 | 2/8 | 13,039.8 | 18.109 | 18.683 | 19.914 | 16.02 |
| `dodecahedron_bevel2_relax_gyro` | 542/900/360/1800 | 4/12 | 15,809.7 | 18.083 | 18.690 | 40.282 | 16.02 |
| `truncatedIcosahedron_truncate50d_ambo_dual` | 272/540/270/1080 | 2/5 | 16,020.9 | 18.050 | 18.684 | 38.654 | 16.00 |
| `icosahedron_kis_gyro` | 272/450/180/900 | 2/14 | 14,427.8 | 15.617 | 16.299 | 31.822 | 16.00 |
| `icosidodecahedron_truncate5d_ambo_dual` | 182/360/180/720 | 4/10 | 14,570.6 | 14.380 | 14.940 | 24.548 | 15.99 |

### Per-pixel figures

Window 2577–2592 records 18,490.8 blended pixels/frame,
1.78× quadrant coverage, at 41.79 cycles/blend.
Inclusive `is_mesh_scan` costs 1316.92 cycles per blended pixel.
Shared-counter parenting can hide `filter_blend` on other draw paths; this
window-local figure is not an invented whole-run blend total.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame 0.46/1.70/32.63 us CPU 3.14%
isr_pack          144.0/frame 6.23/7.12/21.14 us CPU 1.64%
isr_dma_submit    144.0/frame 0.58/0.93/9.29 us CPU 0.22%
```

Times are per-call min/weighted-average/max. Pack averages
7.12 µs/call; submit 0.93 µs/call.
The 600-byte image/black-strobe transfer is asynchronous; its 24 MHz SPI
bound including byte framing is 230 µs, not CPU submit time. ISR share
5.00% leaves approximately 59.38 ms foreground CPU per
62.5 ms window. Render already includes ISRs; peak render needs
1.000× speedup to fit the display budget.

## Summary ranking

1. `is_timeline_step` — 22.615 ms/frame, 36.2% of root time (inclusive).
2. `scan_mesh_raster` — 18.215 ms/frame, 29.2% of root time (inclusive).
3. `is_build_scan` — 7.196 ms/frame, 11.5% of root time (inclusive).

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
`HS_PROFILE_WINDOW=16`; `just profile IslamicStars` runs locked build/flash/capture.
Use the Setup reproduction command for this complete cycle and its flags.
