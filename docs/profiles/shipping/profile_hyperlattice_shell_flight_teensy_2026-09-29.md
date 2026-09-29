# HyperLattice on-device profile — sphere-only shell flight — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time supplemental snapshot of the two authored experimental shell presets with the sphere-only renderer. Replaces the 2026-09-28 shell supplement; its raw evidence remains archived. The canonical cubic roster reports remain unchanged.

Paired report: [global O3](../O3/profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

Portable raw captures, validation, summaries and build provenance are in [the evidence directory](../evidence/hyperlattice_shell_flight_2026-09-29/README.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4; live flywheel and DMA ISRs |
| Image | `profile`; newlib-nano, DMA LEDs, experimental geometry enabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, source `5ebb5655c52ddc654b899935af5c5574035b3c82` |
| Method | Four fixed-preset runs, 30 s each, window 16. Camera and palette move. Exact runtime telemetry excludes frame 1 only; scope trees use complete windows after frame 1. No transition or cycle capture. |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 30 16 "-D HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -D HS_PROFILE_PRESET=<index>"` |

Each instrumented capture attests the default Phantasm image, which passes its size/layout gate.

Preset 5: `FLASH: code:91608, data:153244, headers:9100   free for files:1777664` / `RAM1: variables:315008, code:15224, padding:17544   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

Preset 6: `FLASH: code:91608, data:153244, headers:9100   free for files:1777664` / `RAM1: variables:315008, code:15224, padding:17544   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

The sphere renderer landed in `f14c71a22`; no shell/effect/ray code changed between that commit and the captured source above.

Both authored presets use Shells / 3D, radial mode 0, AA 2, speed 0.025 and spherical shells. Their differing active controls are:

| Index | Cell size | Near fade | Far distance | 3D spin | Shell radius |
|---|--:|--:|--:|--:|--:|
| 5 | 0.78625 | 2 | 10.736 | 0.015 | 0.1 |
| 6 | 0.4645 | 0.6 | 5.836 | 0.003 | 0.15 |

## Frame cadence

Per-frame runtime statistics include the final partial window. Setup frame 1 is excluded from averages, peaks, spill counts and denominators.

| Index | Preset | Mean render ms | Peak render ms | Over 62.5 ms | Over 125 ms | Runtime frames | Captured local |
|---|---|--:|--:|--:|--:|---|---|
| 5 | Shell Flight | 121.418 | 🔴 123.553 | 🔴 228/228 (100.00%) | 0/228 | 2–229 | 2026-09-29 08:03 |
| 6 | Shell Close Flight | 120.839 | 🔴 125.186 | 🔴 228/228 (100.00%) | 1/228 | 2–229 | 2026-09-29 08:06 |

One display window is 62.5 ms at 480 RPM; frames exceeding it miss a flip. The quadrant is 144×72 pixels with a shader sampling margin of 146×73 = 10,658 rays. `canvas_buffer_wait` is display-alignment idle.

### Comparison with the prior matching presets

Same COM4 hardware, authored parameter tuples, 30-second moving-camera method and startup exclusion. The prior 2026-09-28 captures used the general shell renderer. Faster rendering advances more camera frames in the same wall time, so the sampled trajectories differ in length. These are not frame-matched traces; ratios summarize captured means. [Prior evidence](../evidence/hyperlattice_shell_flight_2026-09-28/README.md).

| Preset | Prior mean / peak ms | Sphere mean / peak ms | Capture mean ratio | Observed fps before → after |
|---|--:|--:|--:|--:|
| Shell Flight | 139.052 / 144.076 | 121.418 / 123.553 | 1.145× | 5.34 → 8.01 |
| Shell Close Flight | 135.646 / 139.113 | 120.839 / 125.186 | 1.123× | 5.34 → 7.99 |

## Phase-by-phase readout

Each preset is held independently while camera and palette motion continue. Scope trees use each capture’s worst complete window after setup. Leaves show calls/frame and time/call.

### Shell Flight (frames 81–96)

```text
frame                      124.800 ms  74.880 Mcyc 100.0%
  pov_preserve_half          0.138 ms   0.083 Mcyc   0.1% x1 138us/call
  hl_shader_draw           120.954 ms  72.573 Mcyc  96.9% x1 120954us/call
  hl_timeline_step           0.007 ms   0.004 Mcyc   0.0% x1 7us/call
  canvas_clear               0.094 ms   0.056 Mcyc   0.1% x1 94us/call
  canvas_buffer_wait         1.850 ms   1.110 Mcyc   1.5% x1 1850us/call
```

Wall min/avg/max: 122.965/124.799/125.044 ms. Worst-window render 122.950 ms/frame, selected from 13 complete post-startup windows.

Exactness: root 1198078645 cycles / 600 cycles/µs versus measured wall sum 1996799 µs differs by 0.63 ppm.

Startup setup render: 226.200 ms, excluded. Mean live wall time 124.918 ms corresponds to approximately 8.01 fps.

### Shell Close Flight (frames 193–208)

```text
frame                      124.935 ms  74.961 Mcyc 100.0%
  pov_preserve_half          0.143 ms   0.086 Mcyc   0.1% x1 143us/call
  hl_shader_draw           120.485 ms  72.291 Mcyc  96.4% x1 120485us/call
  hl_timeline_step           0.007 ms   0.004 Mcyc   0.0% x1 7us/call
  canvas_clear               0.093 ms   0.056 Mcyc   0.1% x1 93us/call
  canvas_buffer_wait         2.425 ms   1.455 Mcyc   1.9% x1 2425us/call
```

Wall min/avg/max: 124.747/124.934/125.096 ms. Worst-window render 122.510 ms/frame, selected from 13 complete post-startup windows.

Exactness: root 1199374611 cycles / 600 cycles/µs versus measured wall sum 1998959 µs differs by 0.66 ppm.

Startup setup render: 229.531 ms, excluded. Mean live wall time 125.193 ms corresponds to approximately 7.99 fps.

### Per-pixel figures

The shader writes premultiplied pixels directly; no `filter_blend` counter is used. Cycles per evaluated ray use the selected post-startup scope windows.

| Preset | Shader cycles/ray |
|---|--:|
| Shell Flight | 6809.2 |
| Shell Close Flight | 6782.8 |

## Column-ISR / DMA marshaling cost

Same complete windows as the scope trees.

| Preset / ISR | Calls/frame | Min/avg/max µs | CPU share |
|---|--:|--:|--:|
| Shell Flight / isr_wake | 2301.7 | 0.563/1.680/12.975 | 3.09% |
| Shell Flight / isr_pack | 287.7 | 6.231/6.836/9.918 | 1.57% |
| Shell Flight / isr_dma_submit | 287.7 | 0.691/0.953/1.030 | 0.21% |
| Shell Close Flight / isr_wake | 2304.1 | 0.560/1.700/13.676 | 3.13% |
| Shell Close Flight / isr_pack | 288.0 | 6.231/6.995/10.743 | 1.61% |
| Shell Close Flight / isr_dma_submit | 288.0 | 0.743/0.955/1.941 | 0.22% |

Packing marshals LED data on the CPU; submission launches asynchronous DMA. A 600-byte transfer takes approximately 400 µs at 12 MHz. ISR time is already included in render scopes.

ISR CPU shares total 4.87–4.96%, leaving about 59.4 ms foreground time per display window. Peak render needs approximately 1.98× improvement for Shell Flight and 2.00× for Shell Close Flight to fit 62.5 ms.

## Summary ranking

1. `hl_shader_draw` performs the geometry intersections, antialias coverage and layer compositing; it dominates foreground work.
2. `canvas_buffer_wait` aligns completed rendering to a display flip and is excluded from render timing.
3. Canvas clearing, preserving the opposite half, and choreography are the remaining frame overhead.

## Caveats

- The 30-second paths are bounded samples, not exhaustive orientation or parameter searches.
- CYCCNT includes ISR interruptions. Deep per-ray scopes are disabled; there is no per-pixel instrumentation tax.
- Fixed-preset pinning stops automatic preset transitions; it retains per-frame camera and palette motion.
- The sphere-only shell tracer is a specialized cached-flash hot function with per-frame prepared geometry. Global O3 optimizes the entire instrumented image.
- Device serial telemetry does not expose the simulator’s Unfinished Rays parameter, so these timings do not establish complete traversal for every ray.
- The profiled source was clean. ELF hashes, compiler fingerprints, build flags and captured source status are preserved with the evidence.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`, experimental opt-in, and `HS_PROFILE_PRESET=5` or `6`. Use the Setup command; the wrapper locks COM4 through build, attestation, upload and capture.
