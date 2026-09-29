# HyperLattice Triangular experimental on-device profile (2026-09-27)

The fixed Triangular preset exceeds the 62.5 ms frame budget in every measured
live frame. Shipping selective-O3 averages **126.099 ms**, peaking at
**134.620 ms**; global-O3 averages **134.865 ms**, peaking at **141.587 ms**.

This is a supplemental experiment on clean source
`167fb6e02fe6b5e6c9e04db7639ee755ce01951d`, before the Octet implementation.
The user subsequently requested removal of the Triangular, Cosine, and Gyroid
previews; this report preserves the requested historical Triangular measurement.
It does not replace the [shipping](shipping/profile_hyperlattice_teensy_2026-09-29.md)
or [O3](O3/profile_hyperlattice_teensy_2026-09-27.md) two-preset full-cycle
captures. Those use the normal firmware roster; Triangular is opt-in.

| Configuration | Runtime frames | Mean render ms | Peak render ms | Spilled | Mean wall ms |
| --- | --- | ---: | ---: | ---: | ---: |
| Shipping selective-O3 | 2–442 | 126.099 | 134.620 | 441/441 (100%) | 155.273 |
| Global-O3 | 2–367 | 134.865 | 141.587 | 366/366 (100%) | 187.329 |

## Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0 at 600 MHz, COM3; same physical board, sequential runs; flywheel and DMA ISRs live |
| Images | `profile`: -Os base with selective-O3 shader traversal and color/gamut helpers; `profile_o3`: global -O3 -ffast-math reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice at 288×144; single-entry playlist; fixed `experimental-triangular-flight` |
| Method | 70 seconds per image, 16-frame windows, `HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1`, `HS_PROFILE_PRESET=2`; marker `Profile preset: 2/5` is zero based |
| Motion | Preset choreography paused by the fixed-preset harness; continuous camera motion, rotation, and palette evolution remain active |
| Runtime exclusion | Frame 1 excluded from every runtime average, peak, spill numerator and denominator; subsequent frames all retained |
| Scope exclusion | Entire frames 1–16 window excluded from scope trees and ISR calculations; shipping uses 17–432, O3 uses 17–352 |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice <profile-or-profile_o3> 70 16 '-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2'` |

The captures contain one fixed preset, with no transitions or cycle wrap.
Both untouched logs pass the parser's fixed-preset, monotonic-frame,
complete-telemetry, effect-name, and cycle/wall checks. Exact runtime
statistics also include the trailing 10 shipping and 15 O3 frames after the
last complete counter window. No epoch stretch or dwell compression is used.

Image sizes below are the instrumented single-effect images. The wrapper's
paired Phantasm build is the default firmware image; its attestation is not a
full-roster experimental-build size measurement.

| Image | FLASH code | FLASH data | FLASH headers | RAM1 code | RAM1 variables | RAM1 padding | RAM1 free | RAM2 variables / free |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Shipping | 74,176 | 148,788 | 8,460 | 20,296 | 315,008 | 12,472 | 176,512 | 520,064 / 4,224 |
| O3 | 84,616 | 149,016 | 9,056 | 28,504 | 315,008 | 4,264 | 176,512 | 520,064 / 4,224 |
| O3 minus shipping | +10,440 | +228 | +596 | +8,208 | 0 | -8,208 | 0 | 0 / 0 |

Build logs preserve third-party compiler warnings; these captures do not make
a zero-warning build claim. Compiler and image hashes are in each provenance
file. Build logs and environment dumps are JSON-wrapped with their original
UTF-8 content, line endings, and SHA-256, preserving the exact input bytes.

## Shipping selective-O3

[Raw capture](evidence/hyperlattice_triangular_2026-09-27/ship.txt),
[provenance](evidence/hyperlattice_triangular_2026-09-27/ship.provenance),
[build](evidence/hyperlattice_triangular_2026-09-27/ship_build.json),
[summary](evidence/hyperlattice_triangular_2026-09-27/ship_summary.json),
[parser validation](evidence/hyperlattice_triangular_2026-09-27/ship_validate.txt).
Captured 2026-09-27 21:59 local time (raw-log mtime).

Exactness cross-check: frames 305–320, root
1,199,633,865 cycles / 600 MHz versus 1,999,390 μs wall sum:
**0.113 ppm**.

### Frame cadence

Runtime frames **2–442**: render minimum/mean/peak
**120.231/126.099/134.620 ms**,
spilled **441/441 (100%)**. Startup frame 1 renders in
**229.849 ms**, separately excluded as construction/setup work.
Mean wall time is 155.273 ms, or approximately
6.44 rendered frames per second.

A display interval is 62.5 ms. The nominal quadrant is 144×72 = 10,368 pixels;
the one-pixel shader margin evaluates 146×73 = 10,658 samples. The
`canvas_buffer_wait` scope is alignment idle to a display flip.
Shipping alternates between roughly two-interval (8 fps) and three-interval
(5.33 fps) cadence as the moving camera crosses the 125 ms render threshold.
Even the minimum observed render misses the original one-interval deadline.

### Phase-by-phase readout

There is one held-preset regime. The selected least and most costly complete
shader windows show its camera-dependent variation; neither includes startup.

#### Least costly shader window (frames 17–32)

```text
frame                       124.902 ms  74.941 Mcyc 100.0%
  pov_preserve_half           0.138 ms   0.083 Mcyc   0.1%
  hl_shader_draw            118.538 ms  71.123 Mcyc  94.9%
  hl_timeline_step            0.005 ms   0.003 Mcyc   0.0%
  canvas_clear                0.097 ms   0.058 Mcyc   0.1%
  canvas_buffer_wait          4.342 ms   2.605 Mcyc   3.5%
```

Wall minimum/mean/maximum: 124.577/124.902/125.044 ms.
Mean render is 120.561 ms. All listed leaf counters run once per
frame, so their frame costs also give milliseconds per call. The shader
includes geometry traversal, coverage, and color; the capture has no finer
breakdown inside it.

#### Most costly shader window (frames 369–384)

```text
frame                       187.507 ms 112.504 Mcyc 100.0%
  pov_preserve_half           0.135 ms   0.081 Mcyc   0.1%
  hl_shader_draw            128.419 ms  77.051 Mcyc  68.5%
  hl_timeline_step            0.006 ms   0.004 Mcyc   0.0%
  canvas_clear                0.088 ms   0.053 Mcyc   0.0%
  canvas_buffer_wait         57.051 ms  34.231 Mcyc  30.4%
```

Wall minimum/mean/maximum: 180.297/187.507/194.700 ms.
Mean render is 130.456 ms. All listed leaf counters run once per
frame, so their frame costs also give milliseconds per call. The shader
includes geometry traversal, coverage, and color; the capture has no finer
breakdown inside it.

### Per-pixel figures

This shader writes pixels directly and has no `filter_blend` calls. Across
frames 17–432, `hl_shader_draw` averages
124.305 ms/frame, or 6997.8 cycles per
evaluated sample. This ratio comes from frame-scope totals and does not add
per-pixel instrumentation. Candidate crossings, contributing layers, and
unfinished-ray counts were not logged in these captures.

### Column-ISR / DMA marshaling cost

```text
isr_wake        2896.9/f  0.611/1.672/19.370 us  3.08%
isr_pack         362.1/f  6.230/6.824/10.761 us  1.57%
isr_dma_submit   362.1/f  0.648/0.955/3.160 us  0.22%
```

Times are per-call minimum/average/maximum, followed by CPU share. The rates
are calls per rendered frame; longer frames contain more display interrupts.

- Pack costs 6.824 μs/call; submit costs
  0.955 μs/call.
- The 72-LED image plus black strobe is 600 bytes; wire transfer at 12 MHz
  takes approximately 400 μs asynchronously in DMA.
- Total ISR share is 4.87%, or 3.046 ms per 62.5 ms display
  interval, leaving approximately 59.454 ms for foreground work.
  Render measurements already include interrupts. The observed mean/peak need
  2.02×/2.15× reduction to fit the display
  interval; do not subtract interrupt time twice.

### Summary ranking

1. `hl_shader_draw`: 124.305 ms/frame, 79.1% of root cycles; analytic
   plane-crossing enumeration and coverage share this scope with shading.
2. `canvas_buffer_wait`: 30.765 ms/frame,
   19.6% of root cycles;
   display synchronization rather than renderer computation.
3. Preserve/clear and unscoped preparation make up the remaining work.

## Global-O3

[Raw capture](evidence/hyperlattice_triangular_2026-09-27/o3.txt),
[provenance](evidence/hyperlattice_triangular_2026-09-27/o3.provenance),
[build](evidence/hyperlattice_triangular_2026-09-27/o3_build.json),
[summary](evidence/hyperlattice_triangular_2026-09-27/o3_summary.json),
[parser validation](evidence/hyperlattice_triangular_2026-09-27/o3_validate.txt).
Captured 2026-09-27 22:01 local time (raw-log mtime).

Exactness cross-check: frames 81–96, root
1,799,993,471 cycles / 600 MHz versus 2,999,989 μs wall sum:
**0.039 ppm**.

### Frame cadence

Runtime frames **2–367**: render minimum/mean/peak
**126.210/134.865/141.587 ms**,
spilled **366/366 (100%)**. Startup frame 1 renders in
**243.096 ms**, separately excluded as construction/setup work.
Mean wall time is 187.329 ms, or approximately
5.34 rendered frames per second.

A display interval is 62.5 ms. The nominal quadrant is 144×72 = 10,368 pixels;
the one-pixel shader margin evaluates 146×73 = 10,658 samples. The
`canvas_buffer_wait` scope is alignment idle to a display flip.
Every observed O3 render exceeds 125 ms, keeping this run near three-interval
cadence (5.33 fps). Global optimization does not recover the frame budget.

### Phase-by-phase readout

There is one held-preset regime. The selected least and most costly complete
shader windows show its camera-dependent variation; neither includes startup.

#### Least costly shader window (frames 17–32)

```text
frame                       187.403 ms 112.442 Mcyc 100.0%
  pov_preserve_half           0.137 ms   0.082 Mcyc   0.1%
  hl_shader_draw            126.052 ms  75.631 Mcyc  67.3%
  hl_timeline_step            0.011 ms   0.006 Mcyc   0.0%
  canvas_clear                0.086 ms   0.052 Mcyc   0.0%
  canvas_buffer_wait         59.355 ms  35.613 Mcyc  31.7%
```

Wall minimum/mean/maximum: 184.867/187.403/189.954 ms.
Mean render is 128.048 ms. All listed leaf counters run once per
frame, so their frame costs also give milliseconds per call. The shader
includes geometry traversal, coverage, and color; the capture has no finer
breakdown inside it.

#### Most costly shader window (frames 81–96)

```text
frame                       187.499 ms 112.500 Mcyc 100.0%
  pov_preserve_half           0.135 ms   0.081 Mcyc   0.1%
  hl_shader_draw            135.635 ms  81.381 Mcyc  72.3%
  hl_timeline_step            0.013 ms   0.008 Mcyc   0.0%
  canvas_clear                0.087 ms   0.052 Mcyc   0.0%
  canvas_buffer_wait         49.864 ms  29.918 Mcyc  26.6%
```

Wall minimum/mean/maximum: 182.153/187.499/192.740 ms.
Mean render is 137.635 ms. All listed leaf counters run once per
frame, so their frame costs also give milliseconds per call. The shader
includes geometry traversal, coverage, and color; the capture has no finer
breakdown inside it.

### Per-pixel figures

This shader writes pixels directly and has no `filter_blend` calls. Across
frames 17–352, `hl_shader_draw` averages
133.077 ms/frame, or 7491.7 cycles per
evaluated sample. This ratio comes from frame-scope totals and does not add
per-pixel instrumentation. Candidate crossings, contributing layers, and
unfinished-ray counts were not logged in these captures.

### Column-ISR / DMA marshaling cost

```text
isr_wake        3456.5/f  0.391/1.545/19.860 us  2.85%
isr_pack         432.0/f  5.990/6.721/9.558 us  1.55%
isr_dma_submit   432.0/f  0.605/0.937/4.155 us  0.22%
```

Times are per-call minimum/average/maximum, followed by CPU share. The rates
are calls per rendered frame; longer frames contain more display interrupts.

- Pack costs 6.721 μs/call; submit costs
  0.937 μs/call.
- The 72-LED image plus black strobe is 600 bytes; wire transfer at 12 MHz
  takes approximately 400 μs asynchronously in DMA.
- Total ISR share is 4.61%, or 2.883 ms per 62.5 ms display
  interval, leaving approximately 59.617 ms for foreground work.
  Render measurements already include interrupts. The observed mean/peak need
  2.16×/2.27× reduction to fit the display
  interval; do not subtract interrupt time twice.

### Summary ranking

1. `hl_shader_draw`: 133.077 ms/frame, 71.0% of root cycles; analytic
   plane-crossing enumeration and coverage share this scope with shading.
2. `canvas_buffer_wait`: 52.386 ms/frame,
   27.9% of root cycles;
   display synchronization rather than renderer computation.
3. Preserve/clear and unscoped preparation make up the remaining work.

## Global -O3 versus selective -O3

The two 70-second captures reach different camera-frame ranges. Matching
frames **2–367** gives 366 live frames per build:

| Image | Mean render ms | Peak render ms | Spilled |
| --- | ---: | ---: | ---: |
| Shipping selective-O3 | 126.229 | 133.272 | 366/366 (100%) |
| Global-O3 | 134.865 | 141.587 | 366/366 (100%) |

Global-O3 takes **6.84% longer** on this matched sequence; shipping/O3 mean
speedup is **0.936×**. Its image adds **10,440 FLASH code bytes** and
**8,208 ITCM bytes**. These measurements identify a shader-path cost, but do
not isolate which instructions or flash/cache behavior cause the difference.

The earlier simulator preview measured 14.08 ms averaged over 24 full-canvas
host frames, including initialization. It uses different hardware, sampling
extent, and frame range; it is not a device deadline or a matched speedup
baseline. The original normal-roster cubic/hypercube captures remain separate.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- No per-pixel profiling is enabled. Direct writes have no `filter_blend`
  parenting artifact or added per-pixel scope overhead.
- Shipping retains selective-O3 shader traversal and color/gamut helpers.
  The experimental shader entry and framework validation use the landed
  `HS_HOT_FLASH_MEMBER` placement; the global-O3 run changes compiler flags,
  not those placement annotations.
- Both runs use clean source and the authored Triangular preset. Other
  parameter values, geometry, later implementations, and arbitrary camera
  paths are outside this measurement.
- The fixed-preset knob pauses preset changes, not camera and palette motion.
  No cycle coverage is claimed; the normal roster reports keep their full-cycle
  results and ranking.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_PRESET=2`, and
`HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1`; both runs use the locked
`tools/profile_one.sh` reproduce command above. `just profile HyperLattice`
without these flags profiles the normal two-preset firmware configuration.
