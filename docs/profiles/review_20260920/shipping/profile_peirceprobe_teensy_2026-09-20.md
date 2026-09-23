# PeirceProbe on-device profile — Teensy 4.0, segmented mode (2026-09-20, **selective -O3**)

Point-in-time experimental snapshot. This nonshipping diagnostic does not replace any roster report. Raw capture: `C:/work/Holosphere/build/prof/review_20260920/115/candidate_ship.log`; captured 2026-09-20T00:23:26 America/Los_Angeles. See [comparison](../peirce_pole_classification.md) for matched A/B/A results and limitations.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 on COM4 @ 600 MHz; flywheel and DMA ISRs live |
| Image | `profile`; -Os + newlib-nano, selective hot regions in shader/scan and HD107S packing; GCC 15.2.1 |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Experimental PeirceProbe 288×144, single-entry playlist; candidate `bf4fa1223439ab092855f10b41259e1d2cd68cdb` |
| Method | HS_PROFILE cycle counters; 70 seconds, 32-frame windows; fixed preset, no cycling or dwell compression |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE=C:/work/Holosphere-profile-115 HS_PROFILE_OUT=<unique-absolute-log> bash tools/profile_one.sh PeirceProbe profile 70 32 '-D HS_PROFILE_PEIRCE_PROBE'` |

Image size:

```text
Memory Usage on Teensy 4.0:
  FLASH: code:70952, data:151396, headers:9076   free for files:1800192
   RAM1: variables:315104, code:26504, padding:6264   free for local variables:176416
   RAM2: variables:520064  free for malloc/new:4224
```

Exactness cross-check: frames 1057–1088: root 1197513437 cycles / 600 = 1995855.728 microseconds versus measured 1995858 microseconds, **1.14 ppm**. Parser validation passes; ELF SHA-256 `ce122901cf9124bc98f0b24e38b0b42f922760eff96a0e894f72be9b39274b71`.

## Frame cadence

Pass aggregate: exact render mean **49.583 ms/frame**, worst window **51.891 ms/frame** (frames 513–544), peak frame render **87.407 ms**, spilled **1/1088 (0.092%)**. The sole spill is the cold first frame. Frames 33–1088 peak at 52.604 ms with no spills.

At 480 RPM a display window is 62.5 ms. Each frame shades one 72×144 quadrant, 10,368 pixels. Steady rendering stays within the 16 fps tier. `canvas_buffer_wait` is deliberate idle until the next flip; it is excluded from render time.

## Phase-by-phase readout

Phase schedule: one cold frame followed by a continuously evolving, fixed-parameter preset. No preset wrap is applicable.

### Startup window (frames 1–32)

```text
frame                     62.756 ms 37.654 Mcyc 100.0%
  pov_preserve_half       134.09 us  0.080 Mcyc   0.2% x0.97 138 us/call
  fx_shader_draw          43.230 ms 25.938 Mcyc  68.9% x1 43230 us/call
  fx_prepare_frame         3.866 ms  2.320 Mcyc   6.2% x1 3866 us/call
  fx_advance               2.288 ms  1.373 Mcyc   3.6% x1 2288 us/call
  fx_timeline_step         66.16 us  0.040 Mcyc   0.1% x1 66 us/call
  canvas_clear             90.47 us  0.054 Mcyc   0.1% x1 90 us/call
  canvas_buffer_wait      13.067 ms  7.840 Mcyc  20.8% x1 13067 us/call
```

Wall min/avg/max = 48.591/62.756/87.407 ms. Leaf scopes show calls/frame and per-call microseconds; the startup preservation scope skips the first frame. The startup window includes the cold frame, while steady windows reflect changing camera and palette state. Display wait fills the remaining cadence interval.

### Worst steady window (frames 513–544)

```text
frame                     62.435 ms 37.461 Mcyc 100.0%
  pov_preserve_half       138.50 us  0.083 Mcyc   0.2% x1 138 us/call
  fx_shader_draw          45.309 ms 27.186 Mcyc  72.6% x1 45309 us/call
  fx_prepare_frame         4.071 ms  2.443 Mcyc   6.5% x1 4071 us/call
  fx_advance               2.198 ms  1.319 Mcyc   3.5% x1 2198 us/call
  fx_timeline_step         74.94 us  0.045 Mcyc   0.1% x1 75 us/call
  canvas_clear             86.53 us  0.052 Mcyc   0.1% x1 87 us/call
  canvas_buffer_wait      10.544 ms  6.327 Mcyc  16.9% x1 10544 us/call
```

Wall min/avg/max = 60.795/62.435/64.132 ms. Leaf scopes show calls/frame and per-call microseconds; the startup preservation scope skips the first frame. The startup window includes the cold frame, while steady windows reflect changing camera and palette state. Display wait fills the remaining cadence interval.

### Per-pixel figures

Direct shader output covers 10,368 pixels/frame; no `filter_blend` counter is used. Steady `fx_shader_draw` averages 43.078 ms/frame, or 2492.9 cycles/pixel including raster dispatch and interrupts. This is not an isolated projection-call measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake          1159.3/f  0.46/ 1.70/25.00 us 3.12%
  isr_pack         143.8/f  6.27/ 6.99/14.24 us 1.60%
  isr_dma_submit   143.8/f  0.62/ 0.96/11.40 us 0.22%
```

Columns show calls/frame, min/mean/max per call, and CPU share over capture windows. Pack and submit are nested inside wake; do not add their shares to wake. Pack costs more than submit; the latter schedules asynchronous DMA. The 600-byte composite for 72 LEDs takes 200 microseconds at 24 MHz SPI; this is wire time, not CPU occupancy. `isr_wake` includes pack and submit, so its measured 3.12% is the inclusive ISR CPU share, about 1.95 ms per 62.5 ms display interval. Steady frames need no speedup to hold their cadence; ISR time is already included in render counters.

## Summary ranking

1. `fx_shader_draw`: 43.078 ms/frame, 69.0% of wall frame.
2. `fx_prepare_frame`: 3.985 ms/frame, 6.4% of wall frame.
3. `fx_advance`: 2.208 ms/frame, 3.5% of wall frame.

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
