# GridWorms on-device profile — Teensy 4.0, segmented mode (2026-10-31, **selective -O3**)

Raw capture: `C:/work/temp/gridworms-playback-study/qualified.log`, captured
**2026-10-31 21:15 America/Los_Angeles** on COM3.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, COM3, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridWorms, 288×144, default `grid-worms` preset and temporal sampling |
| Source | `19a2a7e52772f1a694f587092248cb7e59068ce5` with the pending resource-identity and test-precision fixes |
| Method | 70-second capture, 32-frame counter windows, 34 complete windows, no epoch reset |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 HS_DEVICE_WAIT=300 HS_PIO_JOBS=4 bash tools/profile_one.sh GridWorms profile 70 32` |

The profile wrapper built and attested both images under the supported board
lock. The source tree had the provider-identity and test-precision fixes pending; the retained artifact
directory contains `source_status.txt` and `source.diff`. The rendering equations
and authored preset were unchanged from the named commit. Existing resource
family values and provider template identities were restored by appending Twist. The one-hour profile
epoch exceeds the capture and attachment slack. No preset or timing overrides
were used.

Profile ELF SHA-256: `0e0997d84ca1948318670cb9853b5f0626648a7048eb3b5344d61b9618b1d9a5`.
Phantasm ELF SHA-256: `056d27067addbc5a1d329d5d153a90e4e18f0d3a624e6796349939976c050e8e`.
Artifacts: `C:/work/temp/gridworms-playback-study/artifacts/qualified_19a2a7e52772_0e0997d84ca1`.

The complete 37-effect Phantasm image passes existing size/layout gates:

| Bytes | 36-effect baseline | 37 effects | Delta |
|---|--:|--:|--:|
| RAM1 code | 195,784 | 196,184 | +400 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH data | 779,156 | 780,740 | +1,584 |

The 196,608-byte ITCM ceiling retains 424 bytes. Stack headroom remains 13,920
bytes above the existing 12,288-byte floor. No memory budget was changed.

Exactness cross-check: frames 897-928 root cycles divided
by 600 MHz match the measured wall sum within **4.698 ppm**.

## Frame cadence

**Steady render peak: 56.505 ms**, frame 922; mean
**41.986 ms**;
**0/1087 spilled frames (0.00%)** and **0/1087 at or above 58 ms**.
The steady margin below 58 ms is **1.495 ms** in this capture.

**Setup frame 1: 60.735 ms**, also below the 62.5 ms display
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
frame                     61.392 ms 36.835 Mcyc 100.0%
  pov_preserve_half        0.138 ms  0.083 Mcyc   0.2%
  fx_shader_draw          39.959 ms 23.975 Mcyc  65.1%
    fx_prepare_column      1.289 ms  0.774 Mcyc   2.1%
  fx_prepare_frame         0.520 ms  0.312 Mcyc   0.8%
  fx_advance               0.945 ms  0.567 Mcyc   1.5%
    palette_cycle_step     0.871 ms  0.522 Mcyc   1.4%
  fx_timeline_step         0.063 ms  0.038 Mcyc   0.1%
  canvas_clear             0.091 ms  0.054 Mcyc   0.1%
  canvas_buffer_wait      19.488 ms 11.693 Mcyc  31.7%
```

Wall min/mean/max: **31.164/61.392/83.659 ms**. Scopes are inclusive averages;
children overlap their parent. The prepared-column scope runs 144 times per frame.

### Peak steady-render window (frames 897-928)

```text
frame                     62.434 ms 37.461 Mcyc 100.0%
  pov_preserve_half        0.147 ms  0.088 Mcyc   0.2%
  fx_shader_draw          41.615 ms 24.969 Mcyc  66.7%
    fx_prepare_column      1.354 ms  0.812 Mcyc   2.2%
  fx_prepare_frame         0.492 ms  0.295 Mcyc   0.8%
  fx_advance               0.940 ms  0.564 Mcyc   1.5%
    palette_cycle_step     0.914 ms  0.549 Mcyc   1.5%
  fx_timeline_step         0.064 ms  0.038 Mcyc   0.1%
  canvas_clear             0.088 ms  0.053 Mcyc   0.1%
  canvas_buffer_wait      18.908 ms 11.345 Mcyc  30.3%
```

Wall min/mean/max: **37.023/62.433/87.837 ms**. Scopes are inclusive averages;
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
isr_wake        1152.38/frame  0.515/1.946/29.930 us  3.59%
isr_pack         144.03/frame  6.990/8.816/11.463 us  2.03%
isr_dma_submit   144.03/frame  0.621/0.938/1.711 us  0.21%
```

`isr_wake` includes packing and submission; those shares must not be summed.
The submission row measures CPU work, not asynchronous wire-transfer duration.
DMA-completion and other interrupts are not separately measured. Render scopes
already include interrupt time, so the measured 1.495 ms steady margin does not
require subtracting the ISR rows again.

## Summary ranking

The following inclusive costs describe frames 897-928.
The worst individual render frame is reported separately above.

1. `fx_shader_draw`: 41.615 ms/frame, 66.65% of the root window.
2. `fx_prepare_column`: 1.354 ms/frame, 2.17% of the root window.
3. `fx_advance`: 0.940 ms/frame, 1.51% of the root window.

README cells: peak 🟢 56.505, spilled 🟢 0/1087 (0.00%).

README setup frame 1: 60.735 ms.

The capture has zero spills and already meets the requested steady 58 ms
target. The earlier admission capture at `556ecad04` measured 56.461 ms
steady peak and 60.872 ms setup with zero spills; `baseline.log` and its
artifacts remain adjacent to this repeat capture. No performance optimization or authored parameter change was needed.
No matched native/WASM performance comparison was taken.

## Caveats

- This finite capture qualifies the default preset on COM3; it is not a bound
  over arbitrary edited parameters or an unbounded animation interval.
- Cycle scopes include interrupts and instrumentation overhead. Per-pixel
  `filter_blend` instrumentation is disabled; no blend-count estimate is made.
- Selective `HS_O3` regions use the shipping configuration. This is not a
  global-O3 reference capture.
- The pending source changes restored provider identities and qualified
  fast-math test precision; the artifact retains their full diff. Profile report and index edits were made after the capture.

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
  [PASS] root cyc/600 vs wall sum within 5 ppm (4.7 ppm, frames 897-928)

VALID: GridWorms
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridWorms`,
`HS_PROFILE_WINDOW=32`. Use the reproduction command above for build, verified
upload and capture under the shared per-board lock.
