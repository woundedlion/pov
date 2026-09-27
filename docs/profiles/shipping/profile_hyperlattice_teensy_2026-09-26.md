# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Replaces the earlier fixed-hypercube report at this path with a full two-preset
cycle, including both dimensional transitions. Historical evidence is retained.
[Raw capture](../evidence/spherical_perspective_2026-09-26/ship_cycle.txt), [provenance](../evidence/spherical_perspective_2026-09-26/ship_cycle.provenance),
[build/size](../evidence/spherical_perspective_2026-09-26/ship_cycle_build.txt), and [summary](../evidence/spherical_perspective_2026-09-26/ship_cycle_summary.json).

## Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0 @ 600 MHz, COM4; flywheel and DMA ISRs live |
| Image | `profile`; -Os base with selective-O3 shader traversal and color/gamut helpers |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144; clean source 8441a47efd4570cc09914b1b43c03c758606f058 |
| Method | 120 s; 16-frame windows; automatic 320-frame holds and 240-frame transitions; no dwell compression or epoch override |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 120 16` |

Image size:

```text
FLASH: code:62,016, data:148,412, headers:8,708, free:1,812,480
RAM1: variables:314,944, code:19,896, padding:12,872, free:176,576
RAM2: variables:520,064, free:4,224
```

Exactness cross-check: frames 1409–1424, root
598,588,220 cycles / 600 MHz versus wall sum 997,647 μs:
**0.033 ppm**. The untouched capture passes the profile parser's
validation, including both preset markers and their wrap.

## Frame cadence

All captured live frames **2–1897**: mean render
**45.975 ms**, peak **55.579 ms**, spilled
**0/1896 (0.00%)**. Mean wall time is 62.443 ms.
Startup frame 1 renders in **77.272 ms** and is excluded from every
runtime statistic and ranking.

The capture includes 9 trailing frame records after its last
complete counter window. Exact frame statistics retain them. Complete-window
frame telemetry covers frames 2–1888 (1887 live frames);
scope trees and ISR summaries use frames **17–1888**,
excluding the entire startup-containing window. Across those scope windows,
`hl_shader_draw` averages 43.988 ms/frame.

One display window is 62.5 ms; the observed peak leaves 6.921 ms.
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

### Held preset (frames 1873–1888)

```text
frame                      62.430 ms 37.458 Mcyc 100.0%
  pov_preserve_half         0.139 ms  0.083 Mcyc   0.2%
  hl_shader_draw           49.613 ms 29.768 Mcyc  79.5%
  hl_timeline_step          0.006 ms  0.004 Mcyc   0.0%
  canvas_clear              0.086 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       10.841 ms  6.504 Mcyc  17.4%
```

Wall min/mean/max: 58.687/62.430/65.836 ms.
The window's mean render is 51.589 ms. Shader evaluation
dominates; the buffer wait is display-alignment idle. All leaf counters in this
tree run once per frame; their listed milliseconds are also cost per call.

### Transition interval (frames 1617–1632)

```text
frame                      62.496 ms 37.498 Mcyc 100.0%
  pov_preserve_half         0.136 ms  0.081 Mcyc   0.2%
  hl_shader_draw           50.643 ms 30.386 Mcyc  81.0%
  hl_timeline_step          0.020 ms  0.012 Mcyc   0.0%
  canvas_clear              0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait        9.871 ms  5.922 Mcyc  15.8%
```

Wall min/mean/max: 59.206/62.496/65.741 ms.
The window's mean render is 52.625 ms. Shader evaluation
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
| 2 | hypercube-flight | 32 | 1873–1888 | 49.613 | 51.589 | 16 |
| 1 | cubic-flight | 37 | 1361–1376 | 42.433 | 44.431 | 16 |

### Per-pixel figures

This is a direct-write shader with no `filter_blend` calls. At the most costly
reported phase window, shader work is 2851.0 cycles per evaluated
sample using the actual 10,658-sample margin-expanded region. This is a
scope-total ratio, not additional per-pixel instrumentation.

## Column-ISR / DMA marshaling cost

Rates are calls/frame; times are per-call minimum/average/maximum.

```text
isr_wake         1152.2/f  0.571/1.675/20.423 us  3.09%
isr_pack          144.0/f  6.230/6.781/9.851 us  1.56%
isr_dma_submit    144.0/f  0.606/0.949/10.113 us  0.22%
```

- Pack and submit cost 6.781 and 0.949 μs per call, respectively.
- The 72-LED image plus black strobe is 600 bytes; at the configured 12 MHz,
  wire transfer is approximately 400 μs and runs asynchronously in DMA.
- Total ISR share is 4.87%, equivalent to 3.044 ms per
  display window and 59.456 ms for foreground work.
  Render counters already include ISR time; do not subtract it twice. No
  speedup is needed for the observed 16 fps cadence.

## Summary ranking

1. `hl_shader_draw`: 43.988 ms/frame, 70.4% of root cycles.
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

Fixed endpoint checks are also retained: [cubic](../evidence/spherical_perspective_2026-09-26/ship_cubic.txt) and
[hypercube](../evidence/spherical_perspective_2026-09-26/ship_hypercube.txt), each captured for 70 s. The cubic run has
zero spills and peak 47.077 ms versus the original baseline's 40.984 ms;
extraction preserves output but does not assert performance equivalence.
The [baseline evidence](../evidence/spherical_periodic_device_2026-09-26/baseline_summary.json)
records both original endpoint checks.
