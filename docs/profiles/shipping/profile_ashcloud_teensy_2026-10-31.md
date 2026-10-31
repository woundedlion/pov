# AshCloud on-device profile â€” Teensy 4.0, segmented mode (2026-10-31, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `build/prof/ashcloud_ship.log`, captured **2026-10-31 03:25**
America/Los_Angeles on COM4. Replaces `profile_ashcloud_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM shader, ripple and noise kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AshCloud 288Ã—144, authored `ash-cloud-v1-copy1` preset, source `4a67716ae975443ceab6e75b1061d537dee4c2d6` |
| Method | `HS_PROFILE`, 32-frame windows, 130 seconds, `HS_PROFILE_EPOCH_REVS=1200`; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_PIO_JOBS=4 bash tools/profile_one.sh AshCloud profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Profile image: `FLASH: code:128864, data:182700, headers:8948`;
`RAM1: variables:313088, code:38408, padding:27128, free:145664`;
`RAM2: variables:520064, free:4224`.

Full Phantasm image: `FLASH: code:1054128, data:779156, headers:8892`;
`RAM1: variables:312736, code:192568, padding:4040, free:14944`;
`RAM2: variables:520064, free:4224`. Existing size/layout gates pass.

| Region | Starting tree `9cfcf4e24` | Final | Delta |
|---|--:|--:|--:|
| RAM1 code | 192,632 | 192,568 | âˆ’64 |
| RAM1 variables | 312,736 | 312,736 | 0 |
| FLASH data | 749,244 | 779,156 | +29,912 |

The ITCM cap remains 196,608 bytes, the stack floor 12,288 bytes and the
persistent-footprint cap 272,384 bytes. No budget was raised. The immutable
matrix and face copies use 4,320 and 24,576 persistent bytes. Compact RGB
midpoint storage first frees 49,152 bytes; ripple validity adds 1,024 bytes.

Exactness cross-check: root cycles/600 MHz match the measured wall sum in
frames 897-928 within **1.029 ppm**. The capture
validator passes effect identity, coverage, chronology, render/wall distinction
and the cycle/wall cross-check.

Profile ELF SHA-256: `0861adf00b9ceb30471cf0840f54c2cbdbe3fd84b35fa30610269b36592f7d77`.
Phantasm ELF SHA-256: `899318fe99300ebbcbc9cac074d1590bf62faecbbf14a241c984fc3e14cdcc82`.
Retained artifacts: `C:/work/temp/ashcloud-20261030/artifacts/ashcloud-integrated_4a67716ae975_0861adf00b9c`. The capture is from a clean
committed tree. Earlier development captures retain their working-tree diffs.

## Frame cadence

**Steady peak: 57.671 ms**, frame 920; mean **50.249 ms**,
**0/2047 spilled (0.00%)**, and **0/2047 at or above 58 ms**.
Frames 2-2048 in complete windows are included; transitions remain
included. Setup frame 1 is **107.306 ms** and spills separately.
The target margin is **0.329 ms** in this finite capture.

README cells: peak 🟢 57.671, spilled 🟢 0/2047 (0.00%).

A display window is 62.5 ms at 480 RPM. AshCloud shades one 10,368-pixel
quadrant at full sample resolution and holds the 16 fps cadence. Wall time
includes display-sync idle; `canvas_buffer_wait` is that idle and is subtracted
from per-frame render telemetry. The `windows` command's all-frame aggregate
includes startup (50.277 ms mean and 1/2048 spilled);
the ranked indexes explicitly separate setup. Trailing raw rows without a
complete counter window are excluded from the ranked aggregate.

A prior 130-second pass on COM3 before integration measured 57.762 ms peak,
50.266 ms mean and 0/2,047 steady spills. This COM4 pass validates the integrated
source on the second board. These observations qualify authored defaults;
unsupported controls use canonical fallbacks and have no sub-58-ms guarantee.

## Phase-by-phase readout

The single authored preset runs continuously. Camera, ripple bursts, projection
walk, palettes and hue clocks keep advancing. There are no discrete preset
markers or per-preset coverage claims.

### Startup-containing window (frames 1-32)

```
frame                     63.219 ms 37.932 Mcyc 100.0%
  pov_preserve_half        0.137 ms  0.082 Mcyc   0.2%
  fx_shader_draw          43.319 ms 25.991 Mcyc  68.5%
    fx_prepare_column      1.633 ms  0.980 Mcyc   2.6% x144 11.34 us/col
  fx_palette_rgb_bake      0.103 ms  0.062 Mcyc   0.2%
  fx_prepare_frame         0.979 ms  0.587 Mcyc   1.5%
  fx_advance               1.102 ms  0.661 Mcyc   1.7%
    palette_cycle_step     0.992 ms  0.595 Mcyc   1.6%
  fx_timeline_step         0.106 ms  0.063 Mcyc   0.2%
  canvas_clear             0.099 ms  0.059 Mcyc   0.2%
  canvas_buffer_wait      16.903 ms 10.142 Mcyc  26.7%
```

Wall min/mean/max: **39.620/63.219/107.307 ms**. Scopes are inclusive; child
rows overlap their parent. Counter-tree percentages here use the root frame.

Startup includes first-use preparation; its frame-1 cost is reported separately.
It does not describe ongoing cadence.

### Worst steady window (frames 897-928)

```
frame                     62.439 ms 37.463 Mcyc 100.0%
  pov_preserve_half        0.135 ms  0.081 Mcyc   0.2%
  fx_shader_draw          50.308 ms 30.185 Mcyc  80.6%
    fx_prepare_column      1.726 ms  1.036 Mcyc   2.8% x144 11.99 us/col
  fx_palette_rgb_bake      0.486 ms  0.292 Mcyc   0.8%
  fx_prepare_frame         1.143 ms  0.686 Mcyc   1.8%
  fx_advance               1.069 ms  0.642 Mcyc   1.7%
    palette_cycle_step     1.018 ms  0.611 Mcyc   1.6%
  fx_timeline_step         0.106 ms  0.064 Mcyc   0.2%
  canvas_clear             0.097 ms  0.058 Mcyc   0.2%
  canvas_buffer_wait       8.631 ms  5.179 Mcyc  13.8%
```

Wall min/mean/max: **57.536/62.438/67.517 ms**. Scopes are inclusive; child
rows overlap their parent. Counter-tree percentages here use the root frame.

This is the worst window mean, containing the steady peak. Shader evaluation
dominates; frame/column preparation and occasional RGB rebakes are smaller.

### Per-pixel figures

There are 10,368 shaded samples per quadrant. The worst window's shader scope
costs 2911.3 cycles
per pixel, including its column preparation and ISR interruptions. This scan
overwrites premultiplied pixels; there is no `filter_blend` count to substitute
for shaded samples.

## Column-ISR / DMA marshaling cost

For frames 897-928, rates are per rendered frame and
times are per call, min/mean/max:

```
isr_wake        1153.66/frame  0.611/2.028/28.583 us  3.73%
isr_pack         144.19/frame  8.846/9.555/12.381 us  2.20%
isr_dma_submit   144.19/frame  0.720/0.941/5.763 us  0.21%
```

Pack and submit run inside wake. **Wake is inclusive: do not add these shares.**
Wire transfer is asynchronous, 200 microseconds for the 600-byte
72-pixel image-plus-black HD107S transfer at 24 MHz. DMA-completion and other interrupts are not
measured by these ISR scopes. The measured wake share leaves about
60.2 ms of foreground CPU in a 62.5 ms window; render telemetry already includes
interruptions, so its 62.5 ms budget must not subtract wake a second time.

## Summary ranking

A separate 70-second PC-sampler capture of the same optimization on COM3
attributes 57.7% of render samples to the fused shader (including inline curl,
folding, projection and color sampling), 19.5% to ripple transfer, 4.3% to the
lattice and 3.0% to scan dispatch. Separate color finalization is 0.1%, down
from 6.2% in the earlier diagnostic. The sampler costs 153 cycles/service on
average; **its timings are not the performance claim**.

The first port of the supplied JSON measured 78.681 ms observed peak and
71.567 ms mean over 671 steady frames. The final capture is longer and includes
later phases. Retained levers are qualified endpoint/column caches, hot-code
placement, exact lookup ownership and specialized kernels. Geometry row
interpolation was rejected after excessive angular error; the sample grid,
driver cadence and existing fidelity bounds were retained.

The source JSON is byte-identical to the supplied file (SHA-256
`42ed2660dc1148fee7253e4167e09b2982c74b82a26888353fede9c8d7dc641e`).
The additional color-cache errors and qualification controls are documented in
`docs/specs/palette_pattern_modulation.md`. The RAM table copies and compact
storage replayed 2,080 RGB16 frames byte-identically against their original
storage paths; two million fold comparisons were byte-identical. Final physical
fast-math qualification observed RGB16 max 59, premultiplied max 34, alpha
5.37e-5 and direction 3.87e-7 under unchanged bounds. Identical-field scalar
color-tail parity passed 4,179,043 assertions. Full-shader calls can differ in
floating-point rounding; universal bitwise geometry parity is not claimed.

The older global-O3 report uses an older shader and remains historical; it is
not a matched optimization comparison. No new shipping/O3 code delta is given.

Integrated validation passes all 125 native CTests, the 120-frame WASM smoke
across 44 effects at both resolutions, and 259 JavaScript tests with 1,639
assertions (one optional test skipped). The affected uploader suite passes
35 tests. Generated presets match all 19 composed documents. Firmware region
and layout checks, source formatting, lint, include and license checks pass.

## Caveats

- DWT scopes include live ISR time. PC percentages are samples, not stage timers.
- `filter_blend` parenting can obscure overlapping trees; this shader does not
  use it. Nested counters are inclusive and must not be summed.
- No per-pixel timing scopes were added; their overhead would distort the loop.
- Shipping selective `HS_O3` regions are used; this is not a global-O3 image.
- The epoch knob extends capture duration, with no dwell compression or change
  to the work per rendered frame.
- The 58 ms result is a finite default-preset observation, with startup separate.
  It is not a universal bound over all controls or arbitrary animation duration.
- Uploads use the verified pinned TyTools 0.9.9 headless backend through the
  shared-device wrapper. The vendor uploader binaries remain original.

## Harness

`targets/Profile/Profile.ino` supplies cycle/window and ISR logging;
`tools/profile_one.sh` locks, builds, attests, uploads and captures. Run
`just profile AshCloud` for the normal capture; use the explicit 130-second
command above to reproduce this longer validation. Raw timing, PC capture,
build output, envdumps and ELF files remain under `C:/work/temp/ashcloud-20261030`.
