# HyperLattice preset 5 on-device profile (2026-09-28, selective-O3)

This is experimental-octet-wide-flight, preset 5 (internal index 4), held with continuous camera and palette motion. This newly added preset supplements the existing Octet views. The 4D preset is unchanged.

[Shipping report](../shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) and [global-O3 report](../O3/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md).

## Setup

| Item | Value |
|---|---|
| Hardware | Teensy 4.0, 600 MHz, COM4, live flywheel and DMA ISRs |
| Driver | POVSegmented<288,4,480>, segment 0 master |
| Build | profile; -Os base, selective O3 shader and octet geometry |
| Source | 6120d4e5eae0e4e83b09de520a0f1d5b2c46ade6; firmware source committed; any report edits preserved in source.diff evidence |
| Captured | 2026-09-28 00:22, America/Los_Angeles |
| Method | 70 seconds, 16-frame windows, held preset; no choreography compression, preset cycle or epoch crossing |
| Geometry | OCTET, THREE_D; sphere radius 1.0, cell size 3.82825, wire radius 0.015 (relative to cell size) |
| Appearance | Near fade 2.0, far distance 10.666, AA strength 2.0, DEPTH color |
| Motion | Speed 0.12750001, 3D spin 0.010155, 4D spin 0 |

[Raw capture](../evidence/hyperlattice_octet_preset5_2026-09-28/octet_preset5_ship.txt), [provenance](../evidence/hyperlattice_octet_preset5_2026-09-28/octet_preset5_ship.provenance), [validation](../evidence/hyperlattice_octet_preset5_2026-09-28/octet_preset5_ship_validate.txt), [summary](../evidence/hyperlattice_octet_preset5_2026-09-28/octet_preset5_ship_summary.json). The evidence directory preserves SHA-256-wrapped build logs and environment dumps. Full ELF/map artifacts are archived locally in C:/work/Holosphere/build/prof/octet_preset3_20260928.

Instrumented single-effect image sizes:

| Region | Bytes |
|---|---:|
| FLASH code | 75224 |
| FLASH data | 148708 |
| FLASH headers | 8516 |
| FLASH free for files | 1799168 |
| RAM1 variables | 315008 |
| RAM1 code | 20376 |
| RAM1 padding | 12392 |
| RAM1 free for local variables | 176512 |
| RAM2 variables | 520064 |
| RAM2 free for malloc/new | 4224 |

Delta from the previous optimized image in the same build configuration: RAM1 code +0 bytes, RAM1 variables +0 bytes, FLASH data +0 bytes, FLASH code +80 bytes.

The harness also builds and attests the default full-roster Phantasm image with profiling compiled out. Its size/layout gates pass.

Exactness cross-check: frames 129–144, 1,199,967,690 root cycles / 600 MHz versus 1,999,946 µs wall sum: **0.075 ppm** difference.

## Frame cadence

| Runtime frames | Mean render ms | Peak render ms | Spilled/live | Observed fps |
|---|---:|---:|---:|---:|
| 2–550 | 73.299 | 76.425 | 549/549 (100.00%) | 8.01 |

Startup frame 1 renders in 132.607 ms and is excluded. Scope and ISR summaries use complete windows over frames 17–544; runtime statistics retain the trailing per-frame telemetry. The display interval is 62.5 ms; canvas_buffer_wait is alignment idle.

## Phase-by-phase readout

The capture has one held-preset regime. Worst shader window: frames 465–480. Each listed scope runs once per frame.

```text
frame                       125.022 ms  75.013 Mcyc 100.0%
  pov_preserve_half           0.138 ms   0.083 Mcyc   0.1%
  hl_shader_draw             73.516 ms  44.110 Mcyc  58.8%
  hl_timeline_step            0.014 ms   0.008 Mcyc   0.0%
  canvas_clear                0.087 ms   0.052 Mcyc   0.1%
  canvas_buffer_wait         49.482 ms  29.689 Mcyc  39.6%
```

Wall minimum/mean/maximum is 124.793/125.022/125.304 ms. Shader traversal and coverage dominate rendering.

### Per-preset figures

| Preset | Complete windows | Shader ms/frame | Evaluated samples/frame | Cycles/sample |
|---|---:|---:|---:|---:|
| 5: experimental-octet-wide-flight | 33 | 71.311 | 10658 | 4014.5 |

The nominal quadrant is 144×72; the shader margin evaluates 146×73 samples. Pixels are written directly, with no filter_blend calls or per-pixel profiling scopes.

## Column-ISR / DMA marshaling cost

```text
isr_wake        2304.1/f 0.575/1.647/21.348 us  3.04%
isr_pack         288.0/f 6.227/6.803/11.012 us  1.57%
isr_dma_submit   288.0/f 0.622/0.945/ 1.257 us  0.22%
```

Columns show calls/frame, minimum/mean/maximum per-call time, and CPU share. Combined ISR share is 4.82%, leaving approximately 59.486 ms foreground time per display interval. Render measurements already include interrupts. Mean/peak render divided by the 62.5 ms interval is 1.173/1.223.

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
  "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=4"
```
