# HyperLattice on-device profile — authored shell flight — Teensy 4.0, segmented mode (2026-09-28, **-O3**)

Point-in-time supplemental snapshot of the two authored experimental shell presets. The canonical cubic roster reports remain unchanged.

Paired report: [shipping selective O3](../shipping/profile_hyperlattice_shell_flight_teensy_2026-09-28.md).

Portable raw captures, validation, summaries and build provenance are in [the evidence directory](../evidence/hyperlattice_shell_flight_2026-09-28/README.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4; live flywheel and DMA ISRs |
| Image | `profile_o3`; newlib-nano, DMA LEDs, experimental geometry enabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, source `4056013a883d869a00aec38ee0f6c949b664f141` |
| Method | Four fixed-preset runs, 30 s each, window 16. Camera and palette move. Exact runtime telemetry excludes frame 1 only; scope trees use complete windows after frame 1. No transition or cycle capture. |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile_o3 30 16 "-D HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -D HS_PROFILE_PRESET=<index>"` |

Image sizes below are instrumented single-effect images. Each capture also attests the default Phantasm image and its passing size/layout gate.

Preset 5: `FLASH: code:105416, data:153344, headers:8504   free for files:1764352` / `RAM1: variables:315008, code:19960, padding:12808   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

Preset 6: `FLASH: code:105416, data:153344, headers:8504   free for files:1764352` / `RAM1: variables:315008, code:19960, padding:12808   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

The renderer matches landed commit `f4eb56258731e829170072aa983b3a8b9c7a6aeb`; the capture provenance records its isolated source commit above.

Default Phantasm size/layout gate: RAM1 code 167,544 B (+256 B against `bbac9c56c`), RAM1 variables 314,784 B (unchanged), FLASH data 735,504 B (+140 B). Stack free remains 12,896 B.

Both authored presets use Shells / 3D, radial mode 0, AA 2, speed 0.025 and stretch 1. Their differing active controls are:

| Index | Cell size | Near fade | Far distance | 3D spin | Shell radius |
|---|--:|--:|--:|--:|--:|
| 5 | 0.78625 | 2 | 10.736 | 0.015 | 0.1 |
| 6 | 0.4645 | 0.6 | 5.836 | 0.003 | 0.15 |

## Frame cadence

Per-frame runtime statistics include every captured live frame, including the final partial window. Startup frame 1 is excluded from both spill counts and denominators.

| Index | Preset | Mean render ms | Peak render ms | Spilled/live frames | Runtime frames | Captured local |
|---|---|--:|--:|--:|---|---|
| 5 | Shell Flight | 137.408 | 🔴 142.252 | 🔴 153/153 (100.00%) | 2–154 | 2026-09-28 23:24 |
| 6 | Shell Close Flight | 133.910 | 🔴 137.024 | 🔴 153/153 (100.00%) | 2–154 | 2026-09-28 23:27 |

The 480 RPM display window is 62.5 ms (16 fps). Render exceeding this window misses the next flip; wall time includes alignment idle in `canvas_buffer_wait`. One quadrant is 144×72 pixels; its shader sampling margin evaluates 146×73 = 10,658 rays.

## Phase-by-phase readout

Each run holds one preset while its camera and palette continue moving. The following trees show each preset’s worst complete post-startup window. Percentages are shares of frame wall cycles.

### Shell Flight (frames 81–96)

```text
frame                      187.316 ms 112.390 Mcyc 100.0%
  pov_preserve_half          0.134 ms   0.080 Mcyc   0.1% x1 134us/call
  hl_shader_draw           137.848 ms  82.709 Mcyc  73.6% x1 137848us/call
  hl_timeline_step           0.004 ms   0.002 Mcyc   0.0% x1 4us/call
  canvas_clear               0.088 ms   0.053 Mcyc   0.0% x1 88us/call
  canvas_buffer_wait        47.498 ms  28.499 Mcyc  25.4% x1 47498us/call
```

Wall min/avg/max: 181.076/187.316/191.846 ms. Worst-window render: 139.819 ms/frame, across 8 complete post-startup windows.

Exactness: root 1798237588 cycles / 600 cycles/µs versus wall sum 2997064 µs differs by 0.45 ppm.

Startup setup render: 254.841 ms; excluded from runtime statistics. Mean wall time 187.122 ms corresponds to approximately 5.34 frames/s.

### Shell Close Flight (frames 113–128)

```text
frame                      187.454 ms 112.472 Mcyc 100.0%
  pov_preserve_half          0.137 ms   0.082 Mcyc   0.1% x1 137us/call
  hl_shader_draw           132.662 ms  79.597 Mcyc  70.8% x1 132662us/call
  hl_timeline_step           0.004 ms   0.002 Mcyc   0.0% x1 4us/call
  canvas_clear               0.088 ms   0.053 Mcyc   0.0% x1 88us/call
  canvas_buffer_wait        52.809 ms  31.685 Mcyc  28.2% x1 52809us/call
```

Wall min/avg/max: 186.278/187.454/188.667 ms. Worst-window render: 134.645 ms/frame, across 8 complete post-startup windows.

Exactness: root 1799558552 cycles / 600 cycles/µs versus wall sum 2999267 µs differs by 0.92 ppm.

Startup setup render: 254.711 ms; excluded from runtime statistics. Mean wall time 187.116 ms corresponds to approximately 5.34 frames/s.

### Per-pixel figures

The shader writes premultiplied pixels directly; there is no `filter_blend` counter or inferred blend coverage. Shader cycles per evaluated ray, from the windows above:

| Preset | Shader cycles/ray |
|---|--:|
| Shell Flight | 7760.3 |
| Shell Close Flight | 7468.3 |

## Column-ISR / DMA marshaling cost

Readings below use the same complete worst windows as the scope trees.

| Preset / ISR | Calls/frame | Min/avg/max µs | CPU share |
|---|--:|--:|--:|
| Shell Flight / isr_wake | 3453.8 | 0.361/1.551/11.975 | 2.86% |
| Shell Flight / isr_pack | 431.7 | 6.008/6.643/8.830 | 1.53% |
| Shell Flight / isr_dma_submit | 431.7 | 0.760/0.936/1.043 | 0.21% |
| Shell Close Flight / isr_wake | 3456.2 | 0.360/1.570/12.041 | 2.89% |
| Shell Close Flight / isr_pack | 432.0 | 6.015/6.808/9.348 | 1.56% |
| Shell Close Flight / isr_dma_submit | 432.0 | 0.751/0.936/1.016 | 0.21% |

Packing performs LED marshaling on the CPU; submission launches asynchronous DMA. A 600-byte transfer takes approximately 400 µs at 12 MHz. ISR work is already included in the measured render scopes.

The sampled ISR shares total 4.60–4.66%, leaving about 59.6 ms of foreground CPU time per display window. Measured peak render would need a 2.28× improvement for Shell Flight and 2.19× for Shell Close Flight to fit 62.5 ms.

## Summary ranking

1. `hl_shader_draw` performs the geometry intersections, antialias coverage and layer compositing; it dominates foreground work.
2. `canvas_buffer_wait` aligns completed rendering to a display flip and is excluded from render timing.
3. Canvas clearing, preserving the opposite half, and choreography are the remaining frame overhead.

## Caveats

- The 30-second paths are bounded samples, not exhaustive orientation or parameter searches.
- CYCCNT includes ISR interruptions. Deep per-ray scopes are disabled; there is no per-pixel instrumentation tax.
- Fixed-preset pinning stops automatic preset transitions; it retains per-frame camera and palette motion.
- Wire and affine tracing use `HS_O3` regions; shell tracing is a specialized cached-flash hot function with per-frame prepared geometry. Global O3 optimizes the whole instrumented image.
- Device serial telemetry does not expose the simulator’s Unfinished Rays parameter, so these timings do not establish complete traversal for every ray.
- The profiled source was clean. ELF hashes, compiler fingerprints, build flags and captured source status are preserved with the evidence.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`, experimental opt-in, and `HS_PROFILE_PRESET=5` or `6`. Use the Setup command; the wrapper locks COM4 through build, attestation, upload and capture.

The capture wrapper was copied to a local diagnostic helper with serial stdout preserved and absolute paths to its unchanged lock/flash helpers. Build, locking, upload, capture verification and provenance checks were unchanged.

Earlier development captures use different settings and are retained only as [historical evidence](../evidence/hyperlattice_shell_flight_2026-09-28/comparisons/README.md); they are not timing comparisons to these authored defaults.

## Global -O3 vs selective -O3

| Preset | Shipping mean ms | O3 mean ms | Mean speedup | FLASH code Δ | RAM1 code Δ |
|---|--:|--:|--:|--:|--:|
| Shell Flight | 139.052 | 137.408 | 1.012× | +13,344 B | +4,736 B |
| Shell Close Flight | 135.646 | 133.910 | 1.013× | +13,344 B | +4,736 B |
