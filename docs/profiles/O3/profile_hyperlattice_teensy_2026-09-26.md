# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-26, **-O3**)

Global-O3 reference for the [shipping profile](../shipping/profile_hyperlattice_teensy_2026-09-26.md).

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Replaces the earlier fixed-hypercube report at this path with a full two-preset
cycle, including both dimensional transitions. Historical evidence is retained.
[Raw capture](../evidence/spherical_perspective_2026-09-26/o3_cycle.txt), [provenance](../evidence/spherical_perspective_2026-09-26/o3_cycle.provenance),
[build/size](../evidence/spherical_perspective_2026-09-26/o3_cycle_build.txt), and [summary](../evidence/spherical_perspective_2026-09-26/o3_cycle_summary.json).

## Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile_o3`; global -O3 -ffast-math single-effect reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144; clean source 8441a47efd4570cc09914b1b43c03c758606f058 |
| Method | 100 s; 16-frame windows; automatic 320-frame holds and 240-frame transitions; no dwell compression or epoch override |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile_o3 100 16` |

Image size:

```text
FLASH: code:71,752, data:148,672, headers:8,952, free:1,802,240
RAM1: variables:314,976, code:28,136, padding:4,632, free:176,544
RAM2: variables:520,064, free:4,224
```

Exactness cross-check: frames 497–512, root
599,696,438 cycles / 600 MHz versus wall sum 999,494 μs:
**0.063 ppm**. The untouched capture passes the profile parser's
validation, including both preset markers and their wrap.

## Frame cadence

All captured live frames **2–1577**: mean render
**44.883 ms**, peak **54.879 ms**, spilled
**0/1576 (0.00%)**. Mean wall time is 62.429 ms.
Startup frame 1 renders in **76.882 ms** and is excluded from every
runtime statistic and ranking.

The capture includes 9 trailing frame records after its last
complete counter window. Exact frame statistics retain them. Complete-window
frame telemetry covers frames 2–1568 (1567 live frames);
scope trees and ISR summaries use frames **17–1568**,
excluding the entire startup-containing window. Across those scope windows,
`hl_shader_draw` averages 42.889 ms/frame.

One display window is 62.5 ms; the observed peak leaves 7.621 ms.
Both holds and transitions sustain 16 fps without render spills. The nominal
display quadrant is 144×72 = 10,368 pixels; the shader's one-pixel margin
evaluates 146×73 = 10,658 samples. `canvas_buffer_wait` is alignment idle to the
next display flip.

## Phase-by-phase readout

The initial cubic hold is followed by alternating hypercube and cubic holds.
Preset markers at frames 320, 879, 1438 announce
the next target; the complete parameter block switches at the transition's
midpoint. There is no dimensional-rift interpolation. The following trees
select the most expensive complete window in each regime.

### Held preset (frames 593–608)

```text
frame                      62.613 ms 37.568 Mcyc 100.0%
  pov_preserve_half         0.138 ms  0.083 Mcyc   0.2%
  hl_shader_draw           49.268 ms 29.561 Mcyc  78.7%
  hl_timeline_step          0.014 ms  0.008 Mcyc   0.0%
  canvas_clear              0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       11.374 ms  6.825 Mcyc  18.2%
```

Wall min/mean/max: 59.210/62.612/66.744 ms.
The window's mean render is 51.238 ms. Shader evaluation
dominates; the buffer wait is display-alignment idle. All leaf counters in this
tree run once per frame; their listed milliseconds are also cost per call.

### Transition interval (frames 465–480)

```text
frame                      62.403 ms 37.442 Mcyc 100.0%
  pov_preserve_half         0.136 ms  0.081 Mcyc   0.2%
  hl_shader_draw           49.350 ms 29.610 Mcyc  79.1%
  hl_timeline_step          0.020 ms  0.012 Mcyc   0.0%
  canvas_clear              0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       11.098 ms  6.659 Mcyc  17.8%
```

Wall min/mean/max: 58.052/62.402/65.722 ms.
The window's mean render is 51.305 ms. Shader evaluation
dominates; the buffer wait is display-alignment idle. All leaf counters in this
tree run once per frame; their listed milliseconds are also cost per call.

### Per-preset table

The clean-hold filter excludes every marked transition-start frame through
240 following frames, plus the startup window. It conservatively omits the
boundary frames rather than mixing two configurations. The marker sequence
2/2 → 1/2 → 2/2 confirms the two-preset wrap. Rows rank by the maximum
clean-hold shader-window mean, not by the broader per-frame cadence buckets.

| # | Preset | Clean windows | Selected frames | Shader ms/frame | Render ms/frame | fps |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 2 | hypercube-flight | 19 | 593–608 | 49.268 | 51.238 | 16 |
| 1 | cubic-flight | 37 | 1361–1376 | 42.340 | 44.331 | 16 |

### Per-pixel figures

This is a direct-write shader with no `filter_blend` calls. At the most costly
reported phase window, shader work is 2778.2 cycles per evaluated
sample using the actual 10,658-sample margin-expanded region. This is a
scope-total ratio, not additional per-pixel instrumentation.

## Column-ISR / DMA marshaling cost

Rates are calls/frame; times are per-call minimum/average/maximum.

```text
isr_wake         1152.2/f  0.326/1.538/22.206 us  2.84%
isr_pack          144.0/f  5.988/6.603/9.655 us  1.52%
isr_dma_submit    144.0/f  0.628/0.954/10.641 us  0.22%
```

- Pack and submit cost 6.603 and 0.954 μs per call, respectively.
- The 72-LED image plus black strobe is 600 bytes; at the configured 12 MHz,
  wire transfer is approximately 400 μs and runs asynchronously in DMA.
- Total ISR share is 4.58%, equivalent to 2.861 ms per
  display window and 59.639 ms for foreground work.
  Render counters already include ISR time; do not subtract it twice. No
  speedup is needed for the observed 16 fps cadence.

## Summary ranking

1. `hl_shader_draw`: 42.889 ms/frame, 68.7% of root cycles.
2. Palette preparation and other unscoped work occupy the remaining render time.
3. `canvas_buffer_wait` is display-alignment idle, not rendering work.

No matched native/WASM timing comparison was collected. The separate periodic
surface experiment measures different geometry and is not an endpoint speedup
comparison.

## Caveats

- Both authored presets and their automatic transitions are covered; arbitrary
  manual parameter extremes are outside this measured admission range.
- All cycle scopes absorb ISR time because CYCCNT free-runs.
- Deep per-pixel profiling is disabled. Direct writes have no `filter_blend`
  parenting artifact or per-pixel cycle-scope overhead.
- Selective-O3 shader traversal and color/gamut helpers remain active in the
  shipping image; effect preparation and selected cold helpers execute from flash.
- No dwell compression or epoch stretch was used. Source was clean at capture;
  later documentation and experiment-header relocation do not rename its source SHA.
- The paired captures use different physical boards, so small differences
  combine compiler configuration, ISR phase, and board variation.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, and
`HS_PROFILE_WINDOW=16`; use the locked reproduce command above. The default
`just profile HyperLattice` is also a locked run, with its default window size.

Global -O3 versus selective -O3: mean render 45.975 →
44.883 ms (1.024×); peak
55.579 → 54.879 ms. Both configurations have zero
runtime spills. Global -O3 adds +9,736 B of FLASH code and
+8,240 B of ITCM to these single-effect images; that is not a full-roster
admission claim.
