# HyperLattice on-device profile — preset 5 (2026-09-28, global-O3)

This is experimental-octet-wide-flight, preset 5 (internal index 4), held with continuous camera and palette motion. This newly added preset supplements the existing Octet views. The 4D preset is unchanged.

[Shipping report](../shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) and [global-O3 report](../O3/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md).

## Setup

| Item | Value |
|---|---|
| Hardware | Teensy 4.0, 600 MHz, COM4, live flywheel and DMA ISRs |
| Driver | POVSegmented<288,4,480>, segment 0 master |
| Build | profile_o3; global -O3 -ffast-math single-effect reference |
| Source | 6120d4e5eae0e4e83b09de520a0f1d5b2c46ade6; firmware source committed; any report edits preserved in source.diff evidence |
| Captured | 2026-09-28 00:24, America/Los_Angeles |
| Method | 70 seconds, 16-frame windows, held preset; no choreography compression, preset cycle or epoch crossing |
| Geometry | OCTET, THREE_D; sphere radius 1.0, cell size 3.82825, wire radius 0.015 (relative to cell size) |
| Appearance | Near fade 2.0, far distance 10.666, AA strength 2.0, DEPTH color |
| Motion | Speed 0.12750001, 3D spin 0.010155, 4D spin 0 |

Raw capture (supporting artifact removed), provenance (supporting artifact removed), validation (supporting artifact removed), summary (supporting artifact removed). Build logs, environment dumps and ELF/map artifacts are local capture outputs, no longer retained in the repository.

Instrumented single-effect image sizes:

| Region | Bytes |
|---|---:|
| FLASH code | 85384 |
| FLASH data | 148916 |
| FLASH headers | 8388 |
| FLASH free for files | 1788928 |
| RAM1 variables | 315008 |
| RAM1 code | 28520 |
| RAM1 padding | 4248 |
| RAM1 free for local variables | 176512 |
| RAM2 variables | 520064 |
| RAM2 free for malloc/new | 4224 |

Delta from the previous optimized image in the same build configuration: RAM1 code +16 bytes, RAM1 variables +0 bytes, FLASH data +0 bytes, FLASH code +104 bytes.

The harness also builds and attests the default full-roster Phantasm image with profiling compiled out. Its size/layout gates pass.

Exactness cross-check: frames 513–528, 1,197,840,645 root cycles / 600 MHz versus 1,996,401 µs wall sum: **0.038 ppm** difference.

## Frame cadence

README cells: peak 🔴 74.925, spilled 🔴 549/549 (100.00%).

| Runtime frames | Mean render ms | Peak render ms | Spilled/live | Observed fps |
|---|---:|---:|---:|---:|
| 2–550 | 71.635 | 74.925 | 549/549 (100.00%) | 8.01 |

Startup frame 1 renders in 129.560 ms and is excluded. Scope and ISR summaries use complete windows over frames 17–544; runtime statistics retain the trailing per-frame telemetry. The display interval is 62.5 ms; canvas_buffer_wait is alignment idle.

## Phase-by-phase readout

The capture has one held-preset regime. Worst shader window: frames 465–480. Each listed scope runs once per frame.

```text
frame                       125.018 ms  75.011 Mcyc 100.0%
  pov_preserve_half           0.137 ms   0.082 Mcyc   0.1%
  hl_shader_draw             71.857 ms  43.114 Mcyc  57.5%
  hl_timeline_step            0.010 ms   0.006 Mcyc   0.0%
  canvas_clear                0.092 ms   0.055 Mcyc   0.1%
  canvas_buffer_wait         51.128 ms  30.677 Mcyc  40.9%
```

Wall minimum/mean/maximum is 124.762/125.018/125.413 ms. Shader traversal and coverage dominate rendering.

### Per-preset figures

| Preset | Complete windows | Shader ms/frame | Evaluated samples/frame | Cycles/sample |
|---|---:|---:|---:|---:|
| 5: experimental-octet-wide-flight | 33 | 69.643 | 10658 | 3920.6 |

The nominal quadrant is 144×72; the shader margin evaluates 146×73 samples. Pixels are written directly, with no filter_blend calls or per-pixel profiling scopes.

## Column-ISR / DMA marshaling cost

```text
isr_wake        2304.1/f 0.392/1.528/20.353 us  2.82%
isr_pack         288.0/f 5.990/6.747/10.808 us  1.55%
isr_dma_submit   288.0/f 0.595/0.939/ 2.220 us  0.22%
```

Columns show calls/frame, minimum/mean/maximum per-call time, and CPU share. Combined ISR share is 4.59%, leaving approximately 59.632 ms foreground time per display interval. Render measurements already include interrupts. Mean/peak render divided by the 62.5 ms interval is 1.146/1.199.

Pack performs CPU-side LED marshaling; submit launches asynchronous DMA. The 600-byte payload occupies about 230 µs of wire time at the requested 24 MHz (LPSPI framing model).

## Summary ranking

1. hl_shader_draw dominates render time: plane events, strut coverage and compositing.
2. canvas_buffer_wait aligns completed frames to the next display flip.
3. Preserve/clear and timeline preparation account for the remaining work.

The previous shipping 3D settings measured 81.407 ms mean and 87.013 ms peak. Changes in this report are preset workload changes, not further code optimizations; the faster motion also samples a different camera trajectory. No matched host timing is claimed.

## Caveats

- Fixed preset pinning leaves camera/palette animation active; there is no full-cycle or transition-cost claim.
- CYCCNT includes ISR time. Startup and its scope window are excluded.
- Shipping uses selective O3 in the shader and octet geometry, executing from cached flash.
- Results cover the recorded trajectory and settings; later camera positions can cost more.
- The experimental preset remains opt-in; the regular cubic roster ranking is unchanged.

## Harness

The existing targets/Profile/Profile.ino harness runs HyperLattice with HS_PROFILE_WINDOW=16. Native unit_hyper_lattice passes with 67094 assertions, including continuous motion for all three Octet presets. Firmware builds, memory gates, and both capture validators pass.

```sh
export HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM4
bash tools/profile_one.sh HyperLattice profile_o3 70 16 \
  "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=4"
```

## Global-O3 comparison

Matched frames 2–550: shipping mean 73.299 ms, global-O3 mean 71.635 ms, ratio 1.023×. O3 minus shipping: FLASH code +10160 bytes; RAM1 code +8144 bytes.
