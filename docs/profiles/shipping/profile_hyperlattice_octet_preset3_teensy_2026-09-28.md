# HyperLattice preset 3 on-device profile (2026-09-28, selective-O3)

This is experimental-octet-flight, preset 3 (internal index 2), held with continuous camera and palette motion. It supersedes the 3D settings in the [optimization snapshot](../hyperlattice_octet_optimization_2026-09-27.md); that snapshot remains historical evidence for its optimization comparisons. This pair was captured before the wide-flight preset was added; source provenance identifies each image. The 4D preset is unchanged.

[Shipping report](../shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) and [global-O3 report](../O3/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md).

## Setup

| Item | Value |
|---|---|
| Hardware | Teensy 4.0, 600 MHz, COM4, live flywheel and DMA ISRs |
| Driver | POVSegmented<288,4,480>, segment 0 master |
| Build | profile; -Os base, selective O3 shader and octet geometry |
| Source | b60fe925b98ef6fecb009ab81cd854e6374ea52a; firmware source committed; any report edits preserved in source.diff evidence |
| Captured | 2026-09-28 00:17, America/Los_Angeles |
| Method | 70 seconds, 16-frame windows, held preset; no choreography compression, preset cycle or epoch crossing |
| Geometry | OCTET, THREE_D; sphere radius 0, cell size 1.74175, wire radius 0.015 (relative to cell size) |
| Appearance | Near fade 2.0, far distance 4.5, AA strength 2.0, DEPTH color |
| Motion | Speed 0.078, 3D spin 0.008265, 4D spin 0 |

[Raw capture](../evidence/hyperlattice_octet_preset3_2026-09-28/octet_preset3_ship.txt), [provenance](../evidence/hyperlattice_octet_preset3_2026-09-28/octet_preset3_ship.provenance), [validation](../evidence/hyperlattice_octet_preset3_2026-09-28/octet_preset3_ship_validate.txt), [summary](../evidence/hyperlattice_octet_preset3_2026-09-28/octet_preset3_ship_summary.json). The evidence directory preserves SHA-256-wrapped build logs and environment dumps. Full ELF/map artifacts are archived locally in C:/work/Holosphere/build/prof/octet_preset3_20260928.

Instrumented single-effect image sizes:

| Region | Bytes |
|---|---:|
| FLASH code | 75176 |
| FLASH data | 148708 |
| FLASH headers | 8564 |
| FLASH free for files | 1799168 |
| RAM1 variables | 315008 |
| RAM1 code | 20376 |
| RAM1 padding | 12392 |
| RAM1 free for local variables | 176512 |
| RAM2 variables | 520064 |
| RAM2 free for malloc/new | 4224 |

Delta from the previous optimized image in the same build configuration: RAM1 code +0 bytes, RAM1 variables +0 bytes, FLASH data +0 bytes, FLASH code +32 bytes.

The harness also builds and attests the default full-roster Phantasm image with profiling compiled out. Its size/layout gates pass.

Exactness cross-check: frames 177–192, 1,199,674,692 root cycles / 600 MHz versus 1,999,458 µs wall sum: **0.090 ppm** difference.

## Frame cadence

| Runtime frames | Mean render ms | Peak render ms | Spilled/live | Observed fps |
|---|---:|---:|---:|---:|
| 2–550 | 69.841 | 76.566 | 549/549 (100.00%) | 8.01 |

Startup frame 1 renders in 129.607 ms and is excluded. Scope and ISR summaries use complete windows over frames 17–544; runtime statistics retain the trailing per-frame telemetry. The display interval is 62.5 ms; canvas_buffer_wait is alignment idle.

## Phase-by-phase readout

The capture has one held-preset regime. Worst shader window: frames 65–80. Each listed scope runs once per frame.

```text
frame                       124.877 ms  74.926 Mcyc 100.0%
  pov_preserve_half           0.138 ms   0.083 Mcyc   0.1%
  hl_shader_draw             70.536 ms  42.322 Mcyc  56.5%
  hl_timeline_step            0.016 ms   0.009 Mcyc   0.0%
  canvas_clear                0.086 ms   0.052 Mcyc   0.1%
  canvas_buffer_wait         52.304 ms  31.382 Mcyc  41.9%
```

Wall minimum/mean/maximum is 124.486/124.877/125.073 ms. Shader traversal and coverage dominate rendering.

### Per-preset figures

| Preset | Complete windows | Shader ms/frame | Evaluated samples/frame | Cycles/sample |
|---|---:|---:|---:|---:|
| 3: experimental-octet-flight | 33 | 67.840 | 10658 | 3819.1 |

The nominal quadrant is 144×72; the shader margin evaluates 146×73 samples. Pixels are written directly, with no filter_blend calls or per-pixel profiling scopes.

## Column-ISR / DMA marshaling cost

```text
isr_wake        2304.1/f 0.573/1.628/20.693 us  3.00%
isr_pack         288.0/f 6.227/6.682/10.825 us  1.54%
isr_dma_submit   288.0/f 0.603/0.945/ 3.263 us  0.22%
```

Columns show calls/frame, minimum/mean/maximum per-call time, and CPU share. Combined ISR share is 4.76%, leaving approximately 59.525 ms foreground time per display interval. Render measurements already include interrupts. Mean/peak render divided by the 62.5 ms interval is 1.117/1.225.

Pack performs CPU-side LED marshaling; submit launches asynchronous DMA. The 600-byte payload occupies about 400 µs of wire time at 12 MHz.

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

## Harness and validation

The existing targets/Profile/Profile.ino harness runs HyperLattice with HS_PROFILE_WINDOW=16. Native unit_hyper_lattice passes with 67094 assertions, including continuous motion for all three Octet presets. Firmware builds, memory gates, and both capture validators pass.

```sh
export HS_PROFILE_TREE="$PWD" HS_TEENSY_PORT=COM4
bash tools/profile_one.sh HyperLattice profile 70 16 \
  "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2"
```
