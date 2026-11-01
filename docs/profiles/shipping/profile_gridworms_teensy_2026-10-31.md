# GridWorms on-device profile — Teensy 4.0, segmented mode (2026-10-31, **selective -O3**)

Raw capture: `C:/work/temp/gridworms-playback-study/baseline.log`, captured
**2026-10-31 20:23 America/Los_Angeles** on COM3.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, COM3, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridWorms, 288×144, default `grid-worms` preset and temporal sampling |
| Source | `556ecad04d21f1ecadccc58add9089c852f31c09` with the 60-second playlist admission diff |
| Method | 70-second capture, 32-frame counter windows, 34 complete windows, no epoch reset |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 HS_DEVICE_WAIT=300 HS_PIO_JOBS=4 bash tools/profile_one.sh GridWorms profile 70 32` |

The profile wrapper built and attested both images under the supported board
lock. The source tree had the playlist addition pending; the retained artifact
directory contains `source_status.txt` and `source.diff`. The effect implementation
and authored preset were unchanged from the named commit. The one-hour profile
epoch exceeds the capture and attachment slack. No preset or timing overrides
were used.

Profile ELF SHA-256: `7577ecab847600ae427b9c0d41ba32b4de007f651623b96556d19f6e8ab123e6`.
Phantasm ELF SHA-256: `48270284f17c6db73591ec6c7897679256f09734c9e730c70493bb49dd52a822`.
Artifacts: `C:/work/temp/gridworms-playback-study/artifacts/baseline_556ecad04d21_7577ecab8476`.

The complete 37-effect Phantasm image passes existing size/layout gates:

| Bytes | 36-effect baseline | 37 effects | Delta |
|---|--:|--:|--:|
| RAM1 code | 195,784 | 196,184 | +400 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH data | 779,156 | 780,740 | +1,584 |

The 196,608-byte ITCM ceiling retains 424 bytes. Stack headroom remains 13,920
bytes above the existing 12,288-byte floor. No memory budget was changed.

Exactness cross-check: frames 897-928 root cycles divided
by 600 MHz match the measured wall sum within **1.423 ppm**.

## Frame cadence

**Steady render peak: 56.461 ms**, frame 922; mean
**41.905 ms**;
**0/1087 spilled frames (0.00%)** and **0/1087 at or above 58 ms**.
The steady margin below 58 ms is **1.539 ms** in this capture.

**Setup frame 1: 60.872 ms**, also below the 62.5 ms display
window. All 1088 frames in complete windows have zero spills. Steady
statistics exclude frame 1 and include frames 2–1088; 8 trailing
raw rows outside complete counter windows are excluded from the ranked table.
The 70-second pass covers more than the 60-second shipping show segment.

At 480 RPM the display window is 62.5 ms and the observed cadence is 16 fps.
`canvas_buffer_wait` measures synchronization idle and is subtracted from the
per-frame render telemetry. Alternating work varies with palette maintenance
and sphere orientation; wall time includes synchronization wait.

## Phase-by-phase readout

The single preset remains active throughout the pass. Startup and the window
containing the largest steady render sample are shown separately.

### Startup (frames 1-32)

```text
frame                     61.395 ms 36.837 Mcyc 100.0%
  pov_preserve_half        0.137 ms  0.082 Mcyc   0.2%
  fx_shader_draw          39.900 ms 23.940 Mcyc  65.0%
    fx_prepare_column      1.247 ms  0.748 Mcyc   2.0%
  fx_prepare_frame         0.518 ms  0.311 Mcyc   0.8%
  fx_advance               0.944 ms  0.566 Mcyc   1.5%
    palette_cycle_step     0.864 ms  0.519 Mcyc   1.4%
  fx_timeline_step         0.066 ms  0.039 Mcyc   0.1%
  canvas_clear             0.091 ms  0.055 Mcyc   0.1%
  canvas_buffer_wait      19.549 ms 11.730 Mcyc  31.8%
```

Wall min/mean/max: **31.093/61.395/83.620 ms**. Scopes are inclusive averages;
children overlap their parent. The prepared-column scope runs 144 times per frame.

### Peak steady-render window (frames 897-928)

```text
frame                     62.438 ms 37.463 Mcyc 100.0%
  pov_preserve_half        0.148 ms  0.089 Mcyc   0.2%
  fx_shader_draw          41.542 ms 24.925 Mcyc  66.5%
    fx_prepare_column      1.320 ms  0.792 Mcyc   2.1%
  fx_prepare_frame         0.492 ms  0.295 Mcyc   0.8%
  fx_advance               0.948 ms  0.569 Mcyc   1.5%
    palette_cycle_step     0.921 ms  0.553 Mcyc   1.5%
  fx_timeline_step         0.064 ms  0.038 Mcyc   0.1%
  canvas_clear             0.088 ms  0.053 Mcyc   0.1%
  canvas_buffer_wait      18.984 ms 11.390 Mcyc  30.4%
```

Wall min/mean/max: **37.006/62.438/87.857 ms**. Scopes are inclusive averages;
children overlap their parent. The prepared-column scope runs 144 times per frame.

The shader scan dominates this window. Its column preparation and the palette
cycle are included in their named parent scopes; their costs must not be added
to those parent totals.

### Per-pixel figures

Each frame shades one 10,368-pixel quadrant. This timing image does not enable
blended-pixel counts, so geometric coverage is not a measured blend count and
no cycles-per-blend figure is inferred.

## Column-ISR / DMA marshaling cost

The peak steady-render counter window carries these separately timestamped ISR
measurements. Calls/frame use its 32-frame interval; CPU shares use
the ISR window's own elapsed time.

```text
isr_wake        1152.44/frame  0.495/1.941/30.191 us  3.58%
isr_pack         144.03/frame  7.056/8.815/11.855 us  2.03%
isr_dma_submit   144.03/frame  0.621/0.941/7.051 us  0.21%
```

`isr_wake` includes packing and submission; those shares must not be summed.
The submission row measures CPU work, not asynchronous wire-transfer duration.
DMA-completion and other interrupts are not separately measured. Render scopes
already include interrupt time, so the measured 1.539 ms steady margin does not
require subtracting the ISR rows again.

## Summary ranking

The following inclusive costs describe frames 897-928.
The worst individual render frame is reported separately above.

1. `fx_shader_draw`: 41.542 ms/frame, 66.53% of the root window.
2. `fx_prepare_column`: 1.320 ms/frame, 2.11% of the root window.
3. `fx_advance`: 0.948 ms/frame, 1.52% of the root window.

README cells: peak 🟢 56.461, spilled 🟢 0/1087 (0.00%).

README setup frame 1: 60.872 ms.

The capture has zero spills and already meets the requested steady 58 ms
target. No performance optimization or authored parameter change was needed.
No matched native/WASM performance comparison was taken.

## Caveats

- This finite capture qualifies the default preset on COM3; it is not a bound
  over arbitrary edited parameters or an unbounded animation interval.
- Cycle scopes include interrupts and instrumentation overhead. Per-pixel
  `filter_blend` instrumentation is disabled; no blend-count estimate is made.
- Selective `HS_O3` regions use the shipping configuration. This is not a
  global-O3 reference capture.
- The pending source change was playlist admission; the artifact retains its
  full diff. Profile report and index edits were made after the capture.

Parser validation:

```text
effect=GridWorms  windows=34
  [PASS] captured >=3 windows (34)
  [PASS] every window names one effect and resolution (GridWorms 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 34 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 34 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (1.4 ppm, frames 897-928)

VALID: GridWorms
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridWorms`,
`HS_PROFILE_WINDOW=32`. Use the reproduction command above for build, verified
upload and capture under the shared per-board lock.
