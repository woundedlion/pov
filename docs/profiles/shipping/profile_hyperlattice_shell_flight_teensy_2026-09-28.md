# HyperLattice on-device profile — authored shell flight — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time supplemental snapshot of the two authored experimental shell presets. The canonical cubic roster reports remain unchanged.

Paired report: [global O3](../O3/profile_hyperlattice_shell_flight_teensy_2026-09-28.md).

Portable raw captures, validation, summaries and build provenance are in [the evidence directory](../evidence/hyperlattice_shell_flight_2026-09-28/README.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4; live flywheel and DMA ISRs |
| Image | `profile`; newlib-nano, DMA LEDs, experimental geometry enabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, source `4056013a883d869a00aec38ee0f6c949b664f141` |
| Method | Four fixed-preset runs, 30 s each, window 16. Camera and palette move. Exact runtime telemetry excludes frame 1 only; scope trees use complete windows after frame 1. No transition or cycle capture. |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 30 16 "-D HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -D HS_PROFILE_PRESET=<index>"` |

Image sizes below are instrumented single-effect images. Each capture also attests the default Phantasm image and its passing size/layout gate.

Preset 5: `FLASH: code:92072, data:153244, headers:8636   free for files:1777664` / `RAM1: variables:315008, code:15224, padding:17544   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

Preset 6: `FLASH: code:92072, data:153244, headers:8636   free for files:1777664` / `RAM1: variables:315008, code:15224, padding:17544   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

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
| 5 | Shell Flight | 139.052 | 🔴 144.076 | 🔴 153/153 (100.00%) | 2–154 | 2026-09-28 23:20 |
| 6 | Shell Close Flight | 135.646 | 🔴 139.113 | 🔴 153/153 (100.00%) | 2–154 | 2026-09-28 23:25 |

The 480 RPM display window is 62.5 ms (16 fps). Render exceeding this window misses the next flip; wall time includes alignment idle in `canvas_buffer_wait`. One quadrant is 144×72 pixels; its shader sampling margin evaluates 146×73 = 10,658 rays.

## Phase-by-phase readout

Each run holds one preset while its camera and palette continue moving. The following trees show each preset’s worst complete post-startup window. Percentages are shares of frame wall cycles.

### Shell Flight (frames 81–96)

```text
frame                      187.298 ms 112.379 Mcyc 100.0%
  pov_preserve_half          0.132 ms   0.079 Mcyc   0.1% x1 132us/call
  hl_shader_draw           139.548 ms  83.729 Mcyc  74.5% x1 139548us/call
  hl_timeline_step           0.006 ms   0.004 Mcyc   0.0% x1 6us/call
  canvas_clear               0.084 ms   0.051 Mcyc   0.0% x1 84us/call
  canvas_buffer_wait        45.763 ms  27.458 Mcyc  24.4% x1 45763us/call
```

Wall min/avg/max: 180.851/187.298/191.914 ms. Worst-window render: 141.536 ms/frame, across 8 complete post-startup windows.

Exactness: root 1798063424 cycles / 600 cycles/µs versus wall sum 2996775 µs differs by 0.88 ppm.

Startup setup render: 257.926 ms; excluded from runtime statistics. Mean wall time 187.114 ms corresponds to approximately 5.34 frames/s.

### Shell Close Flight (frames 113–128)

```text
frame                      187.442 ms 112.465 Mcyc 100.0%
  pov_preserve_half          0.132 ms   0.079 Mcyc   0.1% x1 132us/call
  hl_shader_draw           134.376 ms  80.625 Mcyc  71.7% x1 134376us/call
  hl_timeline_step           0.007 ms   0.004 Mcyc   0.0% x1 7us/call
  canvas_clear               0.085 ms   0.051 Mcyc   0.0% x1 85us/call
  canvas_buffer_wait        51.074 ms  30.644 Mcyc  27.2% x1 51074us/call
```

Wall min/avg/max: 185.952/187.442/188.912 ms. Worst-window render: 136.368 ms/frame, across 8 complete post-startup windows.

Exactness: root 1799443689 cycles / 600 cycles/µs versus wall sum 2999075 µs differs by 0.73 ppm.

Startup setup render: 258.013 ms; excluded from runtime statistics. Mean wall time 187.108 ms corresponds to approximately 5.34 frames/s.

### Per-pixel figures

The shader writes premultiplied pixels directly; there is no `filter_blend` counter or inferred blend coverage. Shader cycles per evaluated ray, from the windows above:

| Preset | Shader cycles/ray |
|---|--:|
| Shell Flight | 7855.9 |
| Shell Close Flight | 7564.8 |

## Column-ISR / DMA marshaling cost

Readings below use the same complete worst windows as the scope trees.

| Preset / ISR | Calls/frame | Min/avg/max µs | CPU share |
|---|--:|--:|--:|
| Shell Flight / isr_wake | 3453.7 | 0.570/1.655/17.241 | 3.05% |
| Shell Flight / isr_pack | 431.7 | 6.231/6.778/9.248 | 1.56% |
| Shell Flight / isr_dma_submit | 431.7 | 0.793/0.955/1.033 | 0.22% |
| Shell Close Flight / isr_wake | 3456.4 | 0.573/1.670/19.095 | 3.07% |
| Shell Close Flight / isr_pack | 432.0 | 6.231/6.918/9.416 | 1.59% |
| Shell Close Flight / isr_dma_submit | 432.0 | 0.773/0.955/1.031 | 0.22% |

Packing performs LED marshaling on the CPU; submission launches asynchronous DMA. A 600-byte transfer takes approximately 400 µs at 12 MHz. ISR work is already included in the measured render scopes.

The sampled ISR shares total 4.83–4.88%, leaving about 59.5 ms of foreground CPU time per display window. Measured peak render would need a 2.31× improvement for Shell Flight and 2.23× for Shell Close Flight to fit 62.5 ms.

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
