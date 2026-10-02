# HyperLattice on-device profile — Octet — Teensy 4.0, segmented mode (2026-09-27, **selective -O3**)

Point-in-time snapshot of the opt-in Octet 3D and 4D presets. These fixed-preset captures supplement the normal HyperLattice roster report; they do not replace its cubic-preset ranking.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0, 600 MHz, COM4; same board throughout, live flywheel and DMA ISRs |
| Image | profile; -Os base with selective-O3 hot shader and octet geometry functions |
| Driver | POVSegmented<288,4,480>, segment 0 master |
| Effect | HyperLattice 288×144, experimental-octet-flight (preset 2) and experimental-octet-4d-slice (preset 3) |
| Method | Separate held-preset captures; 16-frame windows; 70 seconds for 3D, 45 seconds for 4D. Setup frame 1 and its scope window are excluded. Camera and palette motion remain active. No epoch crossing or preset cycle is claimed. |
| Reproduce | HS_PROFILE_TREE=checkout HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile seconds 16 "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=preset" |

Image sizes are instrumented single-effect images. The paired default Phantasm image is separately attested by the harness.

| Preset | FLASH code | FLASH data | Headers | RAM1 code | RAM1 variables | RAM1 free | RAM2 variables/free |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3D | 75144 | 148708 | 8596 | 20376 | 315008 | 176512 | 520064/4224 |
| 4D | 75144 | 148708 | 8596 | 20376 | 315008 | 176512 | 520064/4224 |

## Frame cadence

| Preset | Mean render ms | Peak render ms | Spilled/live | Observed fps |
|---|---:|---:|---:|---:|
| 3D | 81.407 | 87.013 | 549/549 | 8.01 |
| 4D | 177.308 | 188.387 | 231/231 | 5.31 |

A display interval is 62.5 ms. The nominal quadrant is 144×72; the one-pixel shader margin evaluates 146×73 = 10,658 samples. The canvas_buffer_wait scope measures display-alignment idle, not rendering.

## Phase-by-phase readout

There is one held-preset regime per capture. Each tree below is its worst complete shader window after startup; every listed leaf runs once per frame.

### Octet 3D

Raw capture (supporting artifact removed), provenance (supporting artifact removed), validation (supporting artifact removed), summary (supporting artifact removed).

Runtime frames 2–550; startup frame 1 renders in 152.301 ms and is excluded. Scope/ISR summaries use frames 17–544. Trailing live telemetry after the last complete window remains in the runtime statistics.

Exactness cross-check: frames 145–160, 1,199,249,638 root cycles / 600 MHz versus 1,998,750 µs wall sum, **0.302 ppm** difference.

Worst shader window: frames 513–528.

```text
frame                       124.945 ms  74.967 Mcyc 100.0%
  pov_preserve_half           0.140 ms   0.084 Mcyc   0.1%
  hl_shader_draw             84.790 ms  50.874 Mcyc  67.9%
  hl_timeline_step            0.006 ms   0.004 Mcyc   0.0%
  canvas_clear                0.088 ms   0.053 Mcyc   0.1%
  canvas_buffer_wait         38.124 ms  22.874 Mcyc  30.5%
```

Wall minimum/average/maximum: 124.736/124.945/125.048 ms. Geometry traversal and shading dominate. The distinct 4D nearest-edge query costs more than the 3D plane-pair residual calculation.

### Octet 4D

Raw capture (supporting artifact removed), provenance (supporting artifact removed), validation (supporting artifact removed), summary (supporting artifact removed).

Runtime frames 2–232; startup frame 1 renders in 345.485 ms and is excluded. Scope/ISR summaries use frames 17–224. Trailing live telemetry after the last complete window remains in the runtime statistics.

Exactness cross-check: frames 49–64, 1,794,883,333 root cycles / 600 MHz versus 2,991,473 µs wall sum, **0.260 ppm** difference.

Worst shader window: frames 113–128.

```text
frame                       196.025 ms 117.615 Mcyc 100.0%
  pov_preserve_half           0.129 ms   0.077 Mcyc   0.1%
  hl_shader_draw            179.248 ms 107.549 Mcyc  91.4%
  hl_timeline_step            0.018 ms   0.011 Mcyc   0.0%
  canvas_clear                0.089 ms   0.053 Mcyc   0.0%
  canvas_buffer_wait         14.769 ms   8.861 Mcyc   7.5%
```

Wall minimum/average/maximum: 176.862/196.025/249.589 ms. Geometry traversal and shading dominate. The distinct 4D nearest-edge query costs more than the 3D plane-pair residual calculation.

### Per-preset table

Both presets were pinned for their complete captures, so no transition or wrap-to-zero is applicable. Rows rank by average shader cost over complete post-startup windows.

| Preset | Windows | Samples/frame | Shader ms/frame | Render ms/frame | fps |
|---|---:|---:|---:|---:|---:|
| 4D | 13 | 10658 | 175.242 | 177.308 | 5.31 |
| 3D | 33 | 10658 | 79.363 | 81.407 | 8.01 |

### Per-pixel figures

Pixels are written directly; there are no filter_blend calls and no per-pixel profiling overhead.

- 3D: 4467.8 shader cycles per evaluated sample.
- 4D: 9865.4 shader cycles per evaluated sample.

## Column-ISR / DMA marshaling cost

### 3D

```text
isr_wake        2304.4/f 0.575/1.680/20.262 us  3.10%
isr_pack         288.0/f 6.257/7.023/11.002 us  1.62%
isr_dma_submit   288.0/f 0.620/0.945/ 4.690 us  0.22%
```

Per-call minimum/mean/maximum followed by CPU share. Combined ISR share is 4.93%, leaving approximately 59.417 ms foreground time per display interval. Rendering measurements already include interrupts; do not subtract them twice. Mean/peak render require 1.30×/1.39× reduction to fit one display interval.

Pack performs the CPU-side LED marshaling; submit starts asynchronous DMA. The 600-byte LED/strobe payload takes approximately 230 µs at the requested 24 MHz (LPSPI framing model) on the wire.

### 4D

```text
isr_wake        3472.3/f 0.572/1.671/22.708 us  3.08%
isr_pack         434.0/f 6.227/6.812/ 9.853 us  1.57%
isr_dma_submit   434.0/f 0.640/0.945/ 1.257 us  0.22%
```

Per-call minimum/mean/maximum followed by CPU share. Combined ISR share is 4.87%, leaving approximately 59.458 ms foreground time per display interval. Rendering measurements already include interrupts; do not subtract them twice. Mean/peak render require 2.84×/3.01× reduction to fit one display interval.

Pack performs the CPU-side LED marshaling; submit starts asynchronous DMA. The 600-byte LED/strobe payload takes approximately 230 µs at the requested 24 MHz (LPSPI framing model) on the wire.

## Summary ranking

1. hl_shader_draw dominates both presets: plane traversal, nearest-strut distance, coverage and color compositing.
2. canvas_buffer_wait aligns completed images to the next display flip.
3. Preserve/clear, palette stepping and preparation account for the remaining time.

No matched host timing is claimed. The original optimization baseline and every accepted/rejected trial are documented in the optimization ledger (supporting artifact removed).

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
