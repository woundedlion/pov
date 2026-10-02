# HyperLattice on-device profile — Octet — Teensy 4.0, segmented mode (2026-09-27, **-O3**)

Point-in-time snapshot of the opt-in Octet 3D and 4D presets. These fixed-preset captures supplement the normal HyperLattice roster report; they do not replace its cubic-preset ranking.

[Shipping comparison](../shipping/profile_hyperlattice_octet_teensy_2026-09-27.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0, 600 MHz, COM4; same board throughout, live flywheel and DMA ISRs |
| Image | profile_o3; global -O3 -ffast-math reference |
| Driver | POVSegmented<288,4,480>, segment 0 master |
| Effect | HyperLattice 288×144, experimental-octet-flight (preset 2) and experimental-octet-4d-slice (preset 3) |
| Method | Separate held-preset captures; 16-frame windows; 70 seconds for 3D, 45 seconds for 4D. Setup frame 1 and its scope window are excluded. Camera and palette motion remain active. No epoch crossing or preset cycle is claimed. |
| Reproduce | HS_PROFILE_TREE=checkout HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile_o3 seconds 16 "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=preset" |

Image sizes are instrumented single-effect images. The paired default Phantasm image is separately attested by the harness.

| Preset | FLASH code | FLASH data | Headers | RAM1 code | RAM1 variables | RAM1 free | RAM2 variables/free |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3D | 85280 | 148916 | 8492 | 28504 | 315008 | 176512 | 520064/4224 |
| 4D | 85280 | 148916 | 8492 | 28504 | 315008 | 176512 | 520064/4224 |

## Frame cadence

| Preset | Mean render ms | Peak render ms | Spilled/live | Observed fps |
|---|---:|---:|---:|---:|
| 3D | 79.855 | 85.422 | 549/549 | 8.01 |
| 4D | 175.272 | 186.581 | 232/232 | 5.34 |

A display interval is 62.5 ms. The nominal quadrant is 144×72; the one-pixel shader margin evaluates 146×73 = 10,658 samples. The canvas_buffer_wait scope measures display-alignment idle, not rendering.

## Phase-by-phase readout

There is one held-preset regime per capture. Each tree below is its worst complete shader window after startup; every listed leaf runs once per frame.

### Octet 3D

Raw capture (supporting artifact removed), provenance (supporting artifact removed), validation (supporting artifact removed), summary (supporting artifact removed).

Runtime frames 2–550; startup frame 1 renders in 149.489 ms and is excluded. Scope/ISR summaries use frames 17–544. Trailing live telemetry after the last complete window remains in the runtime statistics.

Exactness cross-check: frames 49–64, 1,199,454,576 root cycles / 600 MHz versus 1,999,091 µs wall sum, **0.020 ppm** difference.

Worst shader window: frames 513–528.

```text
frame                       124.949 ms  74.970 Mcyc 100.0%
  pov_preserve_half           0.138 ms   0.083 Mcyc   0.1%
  hl_shader_draw             83.264 ms  49.959 Mcyc  66.6%
  hl_timeline_step            0.010 ms   0.006 Mcyc   0.0%
  canvas_clear                0.089 ms   0.053 Mcyc   0.1%
  canvas_buffer_wait         39.662 ms  23.797 Mcyc  31.7%
```

Wall minimum/average/maximum: 124.728/124.949/125.209 ms. Geometry traversal and shading dominate. The distinct 4D nearest-edge query costs more than the 3D plane-pair residual calculation.

### Octet 4D

Raw capture (supporting artifact removed), provenance (supporting artifact removed), validation (supporting artifact removed), summary (supporting artifact removed).

Runtime frames 2–233; startup frame 1 renders in 342.012 ms and is excluded. Scope/ISR summaries use frames 17–224. Trailing live telemetry after the last complete window remains in the runtime statistics.

Exactness cross-check: frames 193–208, 1,799,019,097 root cycles / 600 MHz versus 2,998,365 µs wall sum, **0.054 ppm** difference.

Worst shader window: frames 113–128.

```text
frame                       188.214 ms 112.929 Mcyc 100.0%
  pov_preserve_half           0.127 ms   0.076 Mcyc   0.1%
  hl_shader_draw            177.142 ms 106.285 Mcyc  94.1%
  hl_timeline_step            0.018 ms   0.011 Mcyc   0.0%
  canvas_clear                0.086 ms   0.051 Mcyc   0.0%
  canvas_buffer_wait          9.120 ms   5.472 Mcyc   4.8%
```

Wall minimum/average/maximum: 176.911/188.214/196.782 ms. Geometry traversal and shading dominate. The distinct 4D nearest-edge query costs more than the 3D plane-pair residual calculation.

### Per-preset table

Both presets were pinned for their complete captures, so no transition or wrap-to-zero is applicable. Rows rank by average shader cost over complete post-startup windows.

| Preset | Windows | Samples/frame | Shader ms/frame | Render ms/frame | fps |
|---|---:|---:|---:|---:|---:|
| 4D | 13 | 10658 | 173.260 | 175.272 | 5.34 |
| 3D | 33 | 10658 | 77.825 | 79.855 | 8.01 |

### Per-pixel figures

Pixels are written directly; there are no filter_blend calls and no per-pixel profiling overhead.

- 3D: 4381.2 shader cycles per evaluated sample.
- 4D: 9753.8 shader cycles per evaluated sample.

## Column-ISR / DMA marshaling cost

### 3D

```text
isr_wake        2304.4/f 0.393/1.565/20.480 us  2.88%
isr_pack         288.0/f 6.030/6.991/11.188 us  1.61%
isr_dma_submit   288.0/f 0.595/0.939/ 2.997 us  0.22%
```

Per-call minimum/mean/maximum followed by CPU share. Combined ISR share is 4.71%, leaving approximately 59.555 ms foreground time per display interval. Rendering measurements already include interrupts; do not subtract them twice. Mean/peak render require 1.28×/1.37× reduction to fit one display interval.

Pack performs the CPU-side LED marshaling; submit starts asynchronous DMA. The 600-byte LED/strobe payload takes approximately 400 µs at 12 MHz on the wire.

### 4D

```text
isr_wake        3455.5/f 0.395/1.587/25.573 us  2.93%
isr_pack         431.9/f 5.988/6.649/ 9.832 us  1.53%
isr_dma_submit   431.9/f 0.650/0.938/ 1.122 us  0.22%
```

Per-call minimum/mean/maximum followed by CPU share. Combined ISR share is 4.67%, leaving approximately 59.579 ms foreground time per display interval. Rendering measurements already include interrupts; do not subtract them twice. Mean/peak render require 2.80×/2.99× reduction to fit one display interval.

Pack performs the CPU-side LED marshaling; submit starts asynchronous DMA. The 600-byte LED/strobe payload takes approximately 400 µs at 12 MHz on the wire.

## Summary ranking

1. hl_shader_draw dominates both presets: plane traversal, nearest-strut distance, coverage and color compositing.
2. canvas_buffer_wait aligns completed images to the next display flip.
3. Preserve/clear, palette stepping and preparation account for the remaining time.

No matched host timing is claimed. The original optimization baseline and every accepted/rejected trial are documented in the optimization ledger (supporting artifact removed).

## Global -O3 vs selective -O3

| Preset | Matched frames | Shipping mean ms | O3 mean ms | Shipping/O3 speedup | FLASH code delta | RAM1 code delta |
|---|---:|---:|---:|---:|---:|---:|
| 3D | 2–550 | 81.407 | 79.855 | 1.019× | +10136 | +8128 |
| 4D | 2–232 | 177.308 | 175.281 | 1.012× | +10136 | +8128 |

## Caveats

- CYCCNT includes ISR time; memory and compiler flags are recorded with each capture.
- The filter_blend parenting artifact is inapplicable to this direct shader. Per-pixel scopes are disabled.
- Shipping hot shader and geometry functions use selective O3 and cached flash; whole-shader ITCM relocation was not used.
- Fixed-preset pinning pauses preset choreography, not camera or palette motion. There is no dwell compression.
- Results cover the authored preset and recorded camera frames. Other slider values and later camera trajectories may cost more.
- Normalized-coordinate reassociation changes floating-point rounding slightly; native oracle, symmetry, boundary, grazing-ray and rendering comparisons constrain the differences.
- Source commits and any working-tree patch are preserved in the harness provenance/artifacts. No unrelated working-tree changes are included.

## Harness

The targets/Profile/Profile.ino harness uses HS_PROFILE_TARGET=HyperLattice and the fixed-preset flags shown above. The default `just profile HyperLattice` command selects the normal roster, so use the explicit experimental flags to reproduce these captures.
