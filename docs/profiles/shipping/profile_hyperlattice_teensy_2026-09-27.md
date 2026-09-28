# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-27, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Replaces the 2026-09-26 standard report with a fresh full-cycle capture of the
same shipping code, before experimental presets are added. Historical evidence
is retained. [Raw capture](../evidence/hyperlattice_2026-09-27/ship_cycle.txt),
[provenance](../evidence/hyperlattice_2026-09-27/ship_cycle.provenance),
[build/size](../evidence/hyperlattice_2026-09-27/ship_cycle_build.txt),
[summary](../evidence/hyperlattice_2026-09-27/ship_cycle_summary.json).

## Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile`; -Os base with selective-O3 shader traversal and color/gamut helpers |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144; single-entry playlist; clean source 6a7d7160e180e55da8cdd85f0d9d0bfee2cdb478 |
| Method | 100 s; 16-frame windows; automatic 320-frame holds and 240-frame transitions; no dwell compression or epoch override; runtime frames 2–1577; scopes 17–1568 |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile 100 16` |

Image size:

```text
FLASH: code:62,016, data:148,412, headers:8,708, free:1,812,480
RAM1: variables:314,944, code:19,896, padding:12,872, free:176,576
RAM2: variables:520,064, free:4,224
```

Exactness cross-check: frames 1393–1408, root
599,538,616 cycles / 600 MHz versus wall sum 999,231 μs:
**0.027 ppm**. The untouched capture passes the profile parser's
validation, including both preset markers and their wrap.

## Frame cadence

All captured live frames **2–1577**: mean render **45.106 ms**,
peak **55.509 ms**, spilled **0/1576 (0.00%)**.
Mean wall time is 62.441 ms. Startup frame 1 renders in
**77.281 ms** and is excluded from every runtime statistic and ranking.

The capture includes 9 trailing frame records after its last complete
counter window. Exact frame statistics retain them. Complete-window frame
telemetry covers frames 2–1568 (1567 live frames);
scope trees and ISR summaries use frames **17–1568**,
excluding the entire startup-containing window. Across those scope windows,
`hl_shader_draw` averages 43.103 ms/frame; its worst selected window is
49.851 ms/frame at frames 465–480.

One display window is 62.5 ms; the observed peak leaves 6.991 ms.
Both holds and transitions sustain 16 fps without render spills. The nominal
display quadrant is 144×72 = 10,368 pixels; the shader's one-pixel margin
evaluates 146×73 = 10,658 samples. `canvas_buffer_wait` is alignment idle to the
next display flip.

### Comparison with the previous profile

The previous source was `8441a47efd4570cc09914b1b43c03c758606f058`.
There is no source diff in `core/`, `effects/`, `hardware/`, `targets/`, or
`platformio.ini` between that capture and this one. This measures repeatability.
Both rows below use the same frame numbers **2–1577** (1576 live frames),
including the same holds and transitions; frame 1 is excluded.

| Capture | Mean render ms | Peak render ms | Spilled |
| --- | ---: | ---: | ---: |
| [2026-09-26](../evidence/spherical_perspective_2026-09-26/ship_cycle.txt) | 45.105 | 55.579 | 0/1576 |
| 2026-09-27 | 45.106 | 55.509 | 0/1576 |

Mean render changed by +0.001 ms;
peak changed by -0.070 ms.
The earlier shipping capture ran 120 s on COM4; this shipping capture runs
100 s on COM3. Frame matching removes the different phase weighting, while
small board and ISR-phase differences remain. Both O3 captures use COM3.
The former full 120 s shipping average of 45.975 ms is not a matched comparison.

## Phase-by-phase readout

The initial cubic hold is followed by alternating hypercube and cubic holds.
Preset markers at frames 320, 879, 1438 announce the next target;
the complete parameter block switches at the transition's midpoint. There is
no dimensional-rift interpolation. Trees select the most expensive complete
window in each regime.

### Held preset (frames 593–608)

```text
frame                      62.642 ms 37.585 Mcyc 100.0%
  pov_preserve_half         0.138 ms  0.083 Mcyc   0.2%
  hl_shader_draw           49.601 ms 29.761 Mcyc  79.2%
  hl_timeline_step          0.006 ms  0.004 Mcyc   0.0%
  canvas_clear              0.086 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       11.084 ms  6.650 Mcyc  17.7%
```

Wall min/mean/max: 59.059/62.642/66.930 ms.
The window's mean render is 51.559 ms. Both regimes are dominated
by shader evaluation; their buffer wait absorbs the remaining display interval.
All leaf counters run once per frame, so their milliseconds also give cost per call.

### Transition interval (frames 465–480)

```text
frame                      62.439 ms 37.464 Mcyc 100.0%
  pov_preserve_half         0.139 ms  0.083 Mcyc   0.2%
  hl_shader_draw           49.851 ms 29.911 Mcyc  79.8%
  hl_timeline_step          0.020 ms  0.012 Mcyc   0.0%
  canvas_clear              0.086 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       10.597 ms  6.358 Mcyc  17.0%
```

Wall min/mean/max: 58.085/62.439/65.765 ms.
The window's mean render is 51.842 ms. Both regimes are dominated
by shader evaluation; their buffer wait absorbs the remaining display interval.
All leaf counters run once per frame, so their milliseconds also give cost per call.

### Per-preset table

The clean-hold filter excludes every marked transition-start frame through
240 following frames, plus the startup window. The marker sequence
2/2 → 1/2 → 2/2 confirms the two-preset wrap. Rows rank by maximum clean-hold
shader-window mean; README buckets include every live transition frame too.

| # | Preset | Clean windows | Selected frames | Shader ms/frame | Render ms/frame | fps |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 2 | hypercube-flight | 19 | 593–608 | 49.601 | 51.559 | 16 |
| 1 | cubic-flight | 37 | 1361–1376 | 42.431 | 44.428 | 16 |

### Per-pixel figures

This direct-write shader has no `filter_blend` calls. At the most costly
reported phase window, shader work is 2806.4 cycles per evaluated
sample using the actual 10,658-sample region. This is a scope-total ratio,
not additional per-pixel instrumentation.

## Column-ISR / DMA marshaling cost

Rates are calls/frame; times are per-call minimum/average/maximum.

```text
isr_wake         1152.2/f  0.571/1.678/22.370 us  3.09%
isr_pack          144.0/f  6.230/6.792/10.043 us  1.57%
isr_dma_submit    144.0/f  0.608/0.949/9.688 us  0.22%
```

- Pack and submit cost 6.792 and 0.949 μs per call, respectively.
- The 72-LED image plus black strobe is 600 bytes; at the configured 12 MHz,
  wire transfer is approximately 400 μs and runs asynchronously in DMA.
- Total ISR share is 4.88%, equivalent to 3.049 ms per
  display window and 59.451 ms for foreground work.
  Render counters already include ISR time; do not subtract it twice.
  No speedup is needed for either observed 16 fps regime.

## Summary ranking

1. `hl_shader_draw`: 43.103 ms/frame, 69.0% of root cycles.
2. Palette preparation and other unscoped work occupy the remaining render time.
3. `canvas_buffer_wait` is display-alignment idle.

No matched native/WASM timing comparison was collected. The separate periodic
surface experiment measures different geometry and is not an endpoint speedup
comparison.

## Caveats

- Both authored presets and their automatic transitions are covered; arbitrary
  manual parameters and subsequently added experimental presets are outside
  this measured range.
- All scopes absorb ISR time because CYCCNT free-runs.
- Deep per-pixel profiling is disabled. Direct writes have no `filter_blend`
  parenting artifact or per-pixel cycle-scope overhead.
- Selective-O3 shader traversal and color/gamut helpers remain active in the
  shipping image; effect preparation and selected cold helpers execute from flash.
- No dwell compression or epoch stretch was used. Source was clean at capture.
  Both current configurations used the same physical board sequentially.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, and
`HS_PROFILE_WINDOW=16`; use the locked reproduce command above. The default
`just profile HyperLattice` is also a locked run, with its default window size.
Captured 2026-09-27 01:02 local time. The source, ELF hashes, and flags are retained
in provenance; environment dumps are stored as JSON preserving their exact
UTF-8 bytes and SHA-256, including original line endings.

Historical fixed-preset checks remain separate from these cycling captures.
Both fixed poses evolve continuously, but their camera histories differ from
the cycling sequence; they cannot be compared directly to the cycle's means.
The retained extraction comparison used 70 s fixed runs and frames 2–1088 in
each endpoint pair, with setup excluded:

| Fixed endpoint | Before extraction mean/peak ms | Extracted mean/peak ms | Spills |
| --- | ---: | ---: | ---: |
| cubic | 37.669/40.984 | 42.806/47.077 | 0/1087 → 0/1087 |
| hypercube | 44.200/50.174 | 49.715/57.797 | 0/1087 → 0/1087 |

[Original baseline evidence](../evidence/spherical_periodic_device_2026-09-26/baseline_summary.json)
records source `79a93c6f8`; the extracted [cubic](../evidence/spherical_perspective_2026-09-26/ship_cubic.txt)
and [hypercube](../evidence/spherical_perspective_2026-09-26/ship_hypercube.txt)
captures record `8441a47efd`. Extraction preserved endpoint output, with a
measurable rendering cost increase. These historical fixed checks were not
recaptured on 2026-09-27.

## Supplemental fixed Triangular experiment

The opt-in `experimental-triangular-flight` preset was measured separately on
clean source `167fb6e02`, with a 70-second fixed-preset capture on COM3.
Runtime frames 2–442 render in **126.099 ms mean**,
**134.620 ms peak**, with **441/441 spills (100%)**.
Startup frame 1 (229.849 ms) is excluded. These results are outside
the normal two-preset roster and do not replace its measurements above.

The [paired experimental report](../hyperlattice_triangular_2026-09-27.md)
contains both image sizes, scope trees, ISR costs, exact frame ranges, matched
comparison, and portable evidence.

## Supplemental Octet 3D single-owner correction

Point-in-time snapshot of the corrected single-owner strut renderer.
This fixed experimental preset is separate from the normal HyperLattice cycle.
Its [shipping](../shipping/profile_hyperlattice_teensy_2026-09-27.md#supplemental-octet-3d-single-owner-correction) and
[global-O3](../O3/profile_hyperlattice_teensy_2026-09-27.md#supplemental-octet-3d-single-owner-correction) captures use the same source.
[Raw capture](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship.txt),
[provenance](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship.provenance),
[summary](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship_summary.json),
[validation](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship_validate.txt).
Captured 2026-09-27 22:50 local time.

### Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0, 600 MHz, COM3, flywheel and DMA ISRs live |
| Image | `profile`, selective -O3; shipping retains HS_O3 shader helpers and HS_HOT_FLASH_MEMBER placement |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | HyperLattice, 288×144, fixed `experimental-octet-flight`, clean source `3c0ad6dcaad361db1dea671933e4e931e91faba0` |
| Method | 70 seconds, 16-frame windows; runtime frames 2–549; setup frame 1 excluded; scope windows 17–544 |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile 70 16 '-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2'` |

```text
teensy_size:   FLASH: code:76184, data:148692, headers:8596   free for files:1798144
teensy_size:    RAM1: variables:315008, code:20376, padding:12392   free for local variables:176512
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Exactness cross-check: frames 433–448, root 1,199,191,741
cycles / 600 MHz versus 1,998,653 μs wall sum: **0.049 ppm**.
Build logs and environment dumps are retained beside the raw capture.

### Frame cadence

Runtime render minimum/mean/peak: **84.503/88.700/96.869 ms**.
Spills: **548/548 (100.0%)**.
Startup frame 1 took 167.999 ms and is excluded.
Mean wall time: 124.847 ms, approximately 8.01 fps.
The display budget is 62.5 ms. Peak render exceeds it by 34.369 ms.
`canvas_buffer_wait` is alignment idle before the next display flip.

### Phase-by-phase readout

The preset is held while its camera and palette continue moving. No preset
cycle or transition coverage is claimed. Worst complete shader window:

#### Fixed Octet 3D (frames 529–544)

```text
frame                       124.868 ms 74.921 Mcyc 100.0%
  pov_preserve_half           0.138 ms  0.083 Mcyc   0.1%
  hl_shader_draw             94.689 ms 56.814 Mcyc  75.8%
  hl_timeline_step            0.010 ms  0.006 Mcyc   0.0%
  canvas_clear                0.089 ms  0.053 Mcyc   0.1%
  canvas_buffer_wait         28.196 ms 16.917 Mcyc  22.6%
```

Wall minimum/mean/maximum: 124.740/124.868/124.943 ms.
Mean render in this window is 96.672 ms. All listed leaf scopes
run once per frame; their frame costs also give milliseconds per call.
The shader includes traversal, coverage and color without a finer breakdown.

#### Per-pixel figures

The nominal quadrant is 144×72 = 10,368 pixels; the shader margin evaluates
146×73 = 10,658 samples. Across the complete runtime windows, the shader
averages 86.633 ms/frame or 4877.1 cycles/sample.
Direct writes have no `filter_blend` calls. Candidate/layer counts were not captured.

### Column-ISR / DMA marshaling cost

```text
isr_wake        2304.4/f 0.575/1.676/24.428 us 3.09%
isr_pack         288.0/f 6.235/6.997/10.885 us 1.61%
isr_dma_submit   288.0/f 0.616/0.944/1.383 us 0.22%
```

Times are per-call minimum/mean/maximum, followed by CPU share.
Pack averages 6.997 μs versus
0.944 μs for submit. The 600-byte LED
image and black strobe take approximately 400 μs asynchronously at 12 MHz.
ISR share totals 4.92%, leaving approximately 59.425 ms
of foreground time per interval. Render already includes interrupts; its
mean/peak need 1.42×/1.55× reduction to fit.

### Summary ranking

1. `hl_shader_draw`: 86.633 ms/frame, 69.4% of root cycles.
2. `canvas_buffer_wait`: 36.277 ms/frame of display synchronization.
3. Preserve, clear and unscoped preparation account for the remainder.

The previous shipping capture averaged 90.406 ms and peaked at
103.733 ms. Matching frame indices 2–549 gives mean render
90.406 ms before and 88.700 ms after: **1.9% less render time**,
or **1.02× speedup**. Camera motion advances per frame;
the equally long captures can reach different frame indices. The previous
capture is source `39225587612b`; the intervening preset/UI and group-storage
changes mean this is a historical comparison, not an isolated rebuild A/B.
Both use the same fixed Octet 3D settings, board, driver, and compiler.


### Caveats

- All scopes include ISR time. No per-pixel profiling overhead is added.
- Direct writes have no `filter_blend` parenting artifact.
- Shipping uses selective-O3 shader traversal and cached-flash placement;
  global-O3 changes compiler flags, not the placement annotations.
- These measurements cover one authored preset and camera-frame range.
  WASM/native correctness tests are not comparable device timings.
- Both images were built from clean source; no engine or instrumentation
  changes were made for these captures. Setup is excluded, with no warmup cut.

### Harness

`targets/Profile/Profile.ino` supplies the existing HS_PROFILE scopes.
Use the locked reproduce command above; `just profile HyperLattice` without
the experimental/fixed-preset flags measures the normal roster.
