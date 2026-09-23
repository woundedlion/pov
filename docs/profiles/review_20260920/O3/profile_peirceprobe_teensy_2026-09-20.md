# PeirceProbe on-device profile — Teensy 4.0, segmented mode (2026-09-20, **-O3**)

Point-in-time experimental snapshot. This nonshipping diagnostic does not replace any roster report. Raw capture: `C:/work/Holosphere/build/prof/review_20260920/115/candidate_o3.log`; captured 2026-09-20T00:25:40 America/Los_Angeles. See [comparison](../peirce_pole_classification.md) for matched A/B/A results and limitations.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 on COM4 @ 600 MHz; flywheel and DMA ISRs live |
| Image | `profile_o3`; global -O3 -ffast-math + newlib-nano; GCC 15.2.1 |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Experimental PeirceProbe 288×144, single-entry playlist; candidate `bf4fa1223439ab092855f10b41259e1d2cd68cdb` |
| Method | HS_PROFILE cycle counters; 70 seconds, 32-frame windows; fixed preset, no cycling or dwell compression |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE=C:/work/Holosphere-profile-115 HS_PROFILE_OUT=<unique-absolute-log> bash tools/profile_one.sh PeirceProbe profile_o3 70 32 '-D HS_PROFILE_PEIRCE_PROBE'` |

Image size:

```text
Memory Usage on Teensy 4.0:
  FLASH: code:86056, data:151508, headers:8196   free for files:1785856
   RAM1: variables:315136, code:38488, padding:27048   free for local variables:143616
   RAM2: variables:520064  free for malloc/new:4224
```

Exactness cross-check: frames 1057–1088: root 1199639514 cycles / 600 = 1999399.190 microseconds versus measured 1999402 microseconds, **1.41 ppm**. Parser validation passes; ELF SHA-256 `0e4b35e07f7e071a562177d9d3e0f721509a3f6b67d43ded4c3f144696c6b450`.

## Frame cadence

Pass aggregate: exact render mean **44.815 ms/frame**, worst window **47.560 ms/frame** (frames 513–544), peak frame render **79.091 ms**, spilled **1/1088 (0.092%)**. The sole spill is the cold first frame. Frames 33–1088 peak at 48.231 ms with no spills.

At 480 RPM a display window is 62.5 ms. Each frame shades one 72×144 quadrant, 10,368 pixels. Steady rendering stays within the 16 fps tier. `canvas_buffer_wait` is deliberate idle until the next flip; it is excluded from render time.

## Phase-by-phase readout

Phase schedule: one cold frame followed by a continuously evolving, fixed-parameter preset. No preset wrap is applicable.

### Startup window (frames 1–32)

```text
frame                     62.392 ms 37.435 Mcyc 100.0%
  pov_preserve_half       134.03 us  0.080 Mcyc   0.2% x0.97 138 us/call
  fx_shader_draw          39.050 ms 23.430 Mcyc  62.6% x1 39050 us/call
  fx_prepare_frame         3.441 ms  2.065 Mcyc   5.5% x1 3441 us/call
  fx_advance               2.289 ms  1.373 Mcyc   3.7% x1 2289 us/call
  fx_timeline_step         78.16 us  0.047 Mcyc   0.1% x1 78 us/call
  canvas_clear             89.44 us  0.054 Mcyc   0.1% x1 89 us/call
  canvas_buffer_wait      17.281 ms 10.369 Mcyc  27.7% x1 17281 us/call
```

Wall min/avg/max = 44.162/62.391/79.091 ms. Leaf scopes show calls/frame and per-call microseconds; the startup preservation scope skips the first frame. The startup window includes the cold frame, while steady windows reflect changing camera and palette state. Display wait fills the remaining cadence interval.

### Worst steady window (frames 513–544)

```text
frame                     62.486 ms 37.491 Mcyc 100.0%
  pov_preserve_half       136.94 us  0.082 Mcyc   0.2% x1 137 us/call
  fx_shader_draw          41.284 ms 24.770 Mcyc  66.1% x1 41284 us/call
  fx_prepare_frame         3.692 ms  2.215 Mcyc   5.9% x1 3692 us/call
  fx_advance               2.246 ms  1.348 Mcyc   3.6% x1 2246 us/call
  fx_timeline_step         88.69 us  0.053 Mcyc   0.1% x1 89 us/call
  canvas_clear             85.50 us  0.051 Mcyc   0.1% x1 86 us/call
  canvas_buffer_wait      14.926 ms  8.955 Mcyc  23.9% x1 14926 us/call
```

Wall min/avg/max = 60.894/62.485/63.994 ms. Leaf scopes show calls/frame and per-call microseconds; the startup preservation scope skips the first frame. The startup window includes the cold frame, while steady windows reflect changing camera and palette state. Display wait fills the remaining cadence interval.

### Per-pixel figures

Direct shader output covers 10,368 pixels/frame; no `filter_blend` counter is used. Steady `fx_shader_draw` averages 38.633 ms/frame, or 2235.7 cycles/pixel including raster dispatch and interrupts. This is not an isolated projection-call measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake          1159.2/f  0.40/ 1.59/28.65 us 2.94%
  isr_pack         143.8/f  0.52/ 6.90/10.14 us 1.58%
  isr_dma_submit   143.8/f  0.58/ 0.94/10.36 us 0.21%
```

Columns show calls/frame, min/mean/max per call, and CPU share over capture windows. Pack and submit are nested inside wake; do not add their shares to wake. Pack costs more than submit; the latter schedules asynchronous DMA. The 600-byte composite for 72 LEDs takes 200 microseconds at 24 MHz SPI; this is wire time, not CPU occupancy. `isr_wake` includes pack and submit, so its measured 2.94% is the inclusive ISR CPU share, about 1.84 ms per 62.5 ms display interval. Steady frames need no speedup to hold their cadence; ISR time is already included in render counters.

## Summary ranking

1. `fx_shader_draw`: 38.633 ms/frame, 61.9% of wall frame.
2. `fx_prepare_frame`: 3.589 ms/frame, 5.7% of wall frame.
3. `fx_advance`: 2.248 ms/frame, 3.6% of wall frame.

No matching native/WASM performance ledger exists for this experimental target. The native checks in the comparison establish correctness only.

## Caveats

- This is a nonshipping workload. The actual Shader effect is simulator-only; no production feature gate was changed.
- Camera rotation → dodecahedral kaleidoscope → PeirceFastSquare → Grid → EdgeFade → triadic palette uses the existing pullback/raster pipeline. The geometry/source constants mirror workbench preset 4; hue speed is clamped to the composed parameter limit.
- A deliberate `0.85 + 0.05 * edge_class` region tint consumes edge metadata so optimization cannot discard the proposed check. It is identical in all arms, but differs from the real workbench rendering.
- All scopes include ISR time. No per-pixel profiling scopes were introduced. `filter_blend` parenting artifacts are inapplicable to this direct shader path.
- Shipping selective-O3 uses existing hot shader/scan regions and HD107S packing; the twin changes global optimization. No dwell-compression knobs were used.
- Every accepted capture records a clean committed source tree and matching firmware provenance. Initial failed Shader/parameter-registration attempts are excluded.
- The device capture proves the check executes on the raster path, not that rare rounded pole-cap inputs occurred. A native 768-sample pole sweep changes 280 edge-class mismatches to zero; candidate unit_projections passes. No claim is made about tiny-cap visual output on this physical raster.

## Harness

`targets/Profile/Profile.ino`, gated experimental the experimental PeirceProbe header restored from the instrumentation patch, `HS_PROFILE_TARGET=PeirceProbe`, `HS_PROFILE_WINDOW=32`. The explicit reproduce command above adds the required experimental flag; a bare `just profile PeirceProbe` is insufficient.

## Global -O3 vs selective -O3

Candidate render means 49.583 / 44.815 ms = 1.106×. Global-O3 FLASH code is +15,104 B and ITCM code +11,984 B versus the candidate selective-O3 image. See the [shipping sibling](../shipping/profile_peirceprobe_teensy_2026-09-20.md).
