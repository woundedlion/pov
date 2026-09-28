# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-27, **-O3**)

[Shipping selective-O3 sibling](../shipping/profile_hyperlattice_teensy_2026-09-27.md).

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Replaces the 2026-09-26 standard report with a fresh full-cycle capture of the
same shipping code, before experimental presets are added. Historical evidence
is retained. [Raw capture](../evidence/hyperlattice_2026-09-27/o3_cycle.txt),
[provenance](../evidence/hyperlattice_2026-09-27/o3_cycle.provenance),
[build/size](../evidence/hyperlattice_2026-09-27/o3_cycle_build.txt),
[summary](../evidence/hyperlattice_2026-09-27/o3_cycle_summary.json).

## Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0 @ 600 MHz, COM3; flywheel and DMA ISRs live |
| Image | `profile_o3`; global -O3 -ffast-math reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144; single-entry playlist; clean source 6a7d7160e180e55da8cdd85f0d9d0bfee2cdb478 |
| Method | 100 s; 16-frame windows; automatic 320-frame holds and 240-frame transitions; no dwell compression or epoch override; runtime frames 2–1577; scopes 17–1568 |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile_o3 100 16` |

Image size:

```text
FLASH: code:71,752, data:148,672, headers:8,952, free:1,802,240
RAM1: variables:314,976, code:28,136, padding:4,632, free:176,544
RAM2: variables:520,064, free:4,224
```

Exactness cross-check: frames 273–288, root
599,188,210 cycles / 600 MHz versus wall sum 998,647 μs:
**0.017 ppm**. The untouched capture passes the profile parser's
validation, including both preset markers and their wrap.

## Frame cadence

All captured live frames **2–1577**: mean render **44.884 ms**,
peak **54.817 ms**, spilled **0/1576 (0.00%)**.
Mean wall time is 62.429 ms. Startup frame 1 renders in
**76.883 ms** and is excluded from every runtime statistic and ranking.

The capture includes 9 trailing frame records after its last complete
counter window. Exact frame statistics retain them. Complete-window frame
telemetry covers frames 2–1568 (1567 live frames);
scope trees and ISR summaries use frames **17–1568**,
excluding the entire startup-containing window. Across those scope windows,
`hl_shader_draw` averages 42.890 ms/frame; its worst selected window is
49.353 ms/frame at frames 465–480.

One display window is 62.5 ms; the observed peak leaves 7.683 ms.
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
| [2026-09-26](../evidence/spherical_perspective_2026-09-26/o3_cycle.txt) | 44.883 | 54.879 | 0/1576 |
| 2026-09-27 | 44.884 | 54.817 | 0/1576 |

Mean render changed by +0.001 ms;
peak changed by -0.062 ms.
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
frame                      62.618 ms 37.571 Mcyc 100.0%
  pov_preserve_half         0.138 ms  0.083 Mcyc   0.2%
  hl_shader_draw           49.265 ms 29.559 Mcyc  78.7%
  hl_timeline_step          0.013 ms  0.008 Mcyc   0.0%
  canvas_clear              0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       11.395 ms  6.837 Mcyc  18.2%
```

Wall min/mean/max: 59.111/62.617/66.796 ms.
The window's mean render is 51.223 ms. Both regimes are dominated
by shader evaluation; their buffer wait absorbs the remaining display interval.
All leaf counters run once per frame, so their milliseconds also give cost per call.

### Transition interval (frames 465–480)

```text
frame                      62.409 ms 37.445 Mcyc 100.0%
  pov_preserve_half         0.137 ms  0.082 Mcyc   0.2%
  hl_shader_draw           49.353 ms 29.612 Mcyc  79.1%
  hl_timeline_step          0.021 ms  0.013 Mcyc   0.0%
  canvas_clear              0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait       11.091 ms  6.655 Mcyc  17.8%
```

Wall min/mean/max: 58.031/62.408/65.724 ms.
The window's mean render is 51.317 ms. Both regimes are dominated
by shader evaluation; their buffer wait absorbs the remaining display interval.
All leaf counters run once per frame, so their milliseconds also give cost per call.

### Per-preset table

The clean-hold filter excludes every marked transition-start frame through
240 following frames, plus the startup window. The marker sequence
2/2 → 1/2 → 2/2 confirms the two-preset wrap. Rows rank by maximum clean-hold
shader-window mean; README buckets include every live transition frame too.

| # | Preset | Clean windows | Selected frames | Shader ms/frame | Render ms/frame | fps |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 2 | hypercube-flight | 19 | 593–608 | 49.265 | 51.223 | 16 |
| 1 | cubic-flight | 37 | 1361–1376 | 42.339 | 44.327 | 16 |

### Per-pixel figures

This direct-write shader has no `filter_blend` calls. At the most costly
reported phase window, shader work is 2778.3 cycles per evaluated
sample using the actual 10,658-sample region. This is a scope-total ratio,
not additional per-pixel instrumentation.

## Column-ISR / DMA marshaling cost

Rates are calls/frame; times are per-call minimum/average/maximum.

```text
isr_wake         1152.2/f  0.326/1.538/22.591 us  2.84%
isr_pack          144.0/f  5.988/6.603/9.671 us  1.52%
isr_dma_submit    144.0/f  0.610/0.954/10.686 us  0.22%
```

- Pack and submit cost 6.603 and 0.954 μs per call, respectively.
- The 72-LED image plus black strobe is 600 bytes; at the configured 12 MHz,
  wire transfer is approximately 400 μs and runs asynchronously in DMA.
- Total ISR share is 4.58%, equivalent to 2.861 ms per
  display window and 59.639 ms for foreground work.
  Render counters already include ISR time; do not subtract it twice.
  No speedup is needed for either observed 16 fps regime.

## Summary ranking

1. `hl_shader_draw`: 42.890 ms/frame, 68.7% of root cycles.
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
Captured 2026-09-27 01:05 local time. The source, ELF hashes, and flags are retained
in provenance; environment dumps are stored as JSON preserving their exact
UTF-8 bytes and SHA-256, including original line endings.

Global -O3 versus selective -O3, both frames 2–1577 on COM3: mean render
45.106 → 44.884 ms (1.005×); peak
55.509 → 54.817 ms. Both have zero runtime spills.
Global -O3 adds +9,736 B of FLASH code and +8,240 B of ITCM to these
single-effect images; that is not a full-roster admission claim.

## Supplemental fixed Triangular experiment

The opt-in `experimental-triangular-flight` preset was measured separately on
clean source `167fb6e02`, with a 70-second fixed-preset capture on COM3.
Runtime frames 2–367 render in **134.865 ms mean**,
**141.587 ms peak**, with **366/366 spills (100%)**.
Startup frame 1 (243.096 ms) is excluded. These results are outside
the normal two-preset roster and do not replace its measurements above.

The [paired experimental report](../hyperlattice_triangular_2026-09-27.md)
contains both image sizes, scope trees, ISR costs, exact frame ranges, matched
comparison, and portable evidence.
