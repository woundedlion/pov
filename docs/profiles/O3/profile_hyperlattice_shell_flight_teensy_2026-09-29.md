# HyperLattice on-device profile — sphere-only shell flight — Teensy 4.0, segmented mode (2026-09-29, **-O3**)

Point-in-time supplemental snapshot of the two authored experimental shell presets with the sphere-only renderer. Replaces the 2026-09-28 shell supplement; its raw evidence remains archived. The canonical cubic roster reports remain unchanged.

Paired report: [shipping selective O3](../shipping/profile_hyperlattice_shell_flight_teensy_2026-09-29.md).

Portable raw captures, validation, summaries and build provenance are in [the evidence directory](../evidence/hyperlattice_shell_flight_2026-09-29/README.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4; live flywheel and DMA ISRs |
| Image | `profile_o3`; newlib-nano, DMA LEDs, experimental geometry enabled |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, source `5ebb5655c52ddc654b899935af5c5574035b3c82` |
| Method | Four fixed-preset runs, 30 s each, window 16. Camera and palette move. Exact runtime telemetry excludes frame 1 only; scope trees use complete windows after frame 1. No transition or cycle capture. |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile_o3 30 16 "-D HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -D HS_PROFILE_PRESET=<index>"` |

Each instrumented capture attests the default Phantasm image, which passes its size/layout gate.

Preset 5: `FLASH: code:104960, data:153344, headers:8960   free for files:1764352` / `RAM1: variables:315008, code:19960, padding:12808   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

Preset 6: `FLASH: code:104960, data:153344, headers:8960   free for files:1764352` / `RAM1: variables:315008, code:19960, padding:12808   free for local variables:176512` / `RAM2: variables:520064  free for malloc/new:4224`.

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
| 5 | Shell Flight | 115.243 | 🔴 117.363 | 🔴 228/228 (100.00%) | 0/228 | 2–229 | 2026-09-29 08:04 |
| 6 | Shell Close Flight | 115.208 | 🔴 121.337 | 🔴 228/228 (100.00%) | 0/228 | 2–229 | 2026-09-29 08:07 |

One display window is 62.5 ms at 480 RPM; frames exceeding it miss a flip. The quadrant is 144×72 pixels with a shader sampling margin of 146×73 = 10,658 rays. `canvas_buffer_wait` is display-alignment idle.

### Comparison with the prior matching presets

Same COM4 hardware, authored parameter tuples, 30-second moving-camera method and startup exclusion. The prior 2026-09-28 captures used the general shell renderer. Faster rendering advances more camera frames in the same wall time, so the sampled trajectories differ in length. These are not frame-matched traces; ratios summarize captured means. [Prior evidence](../evidence/hyperlattice_shell_flight_2026-09-28/README.md).

| Preset | Prior mean / peak ms | Sphere mean / peak ms | Capture mean ratio | Observed fps before → after |
|---|--:|--:|--:|--:|
| Shell Flight | 137.408 / 142.252 | 115.243 / 117.363 | 1.192× | 5.34 → 8.01 |
| Shell Close Flight | 133.910 / 137.024 | 115.208 / 121.337 | 1.162× | 5.34 → 8.01 |

## Phase-by-phase readout

Each preset is held independently while camera and palette motion continue. Scope trees use each capture’s worst complete window after setup. Leaves show calls/frame and time/call.

### Shell Flight (frames 81–96)

```text
frame                      124.822 ms  74.893 Mcyc 100.0%
  pov_preserve_half          0.139 ms   0.083 Mcyc   0.1% x1 139us/call
  hl_shader_draw           114.432 ms  68.659 Mcyc  91.7% x1 114432us/call
  hl_timeline_step           0.003 ms   0.002 Mcyc   0.0% x1 3us/call
  canvas_clear               0.091 ms   0.055 Mcyc   0.1% x1 91us/call
  canvas_buffer_wait         8.261 ms   4.957 Mcyc   6.6% x1 8261us/call
```

Wall min/avg/max: 123.229/124.821/125.052 ms. Worst-window render 116.560 ms/frame, selected from 13 complete post-startup windows.

Exactness: root 1198287732 cycles / 600 cycles/µs versus measured wall sum 1997148 µs differs by 0.89 ppm.

Startup setup render: 214.322 ms, excluded. Mean live wall time 124.908 ms corresponds to approximately 8.01 fps.

### Shell Close Flight (frames 161–176)

```text
frame                      125.156 ms  75.093 Mcyc 100.0%
  pov_preserve_half          0.140 ms   0.084 Mcyc   0.1% x1 140us/call
  hl_shader_draw           116.414 ms  69.849 Mcyc  93.0% x1 116414us/call
  hl_timeline_step           0.004 ms   0.002 Mcyc   0.0% x1 4us/call
  canvas_clear               0.092 ms   0.055 Mcyc   0.1% x1 92us/call
  canvas_buffer_wait         6.621 ms   3.973 Mcyc   5.3% x1 6621us/call
```

Wall min/avg/max: 124.105/125.155/125.810 ms. Worst-window render 118.534 ms/frame, selected from 13 complete post-startup windows.

Exactness: root 1201492947 cycles / 600 cycles/µs versus measured wall sum 2002490 µs differs by 0.88 ppm.

Startup setup render: 218.008 ms, excluded. Mean live wall time 124.909 ms corresponds to approximately 8.01 fps.

### Per-pixel figures

The shader writes premultiplied pixels directly; no `filter_blend` counter is used. Cycles per evaluated ray use the selected post-startup scope windows.

| Preset | Shader cycles/ray |
|---|--:|
| Shell Flight | 6442.0 |
| Shell Close Flight | 6553.6 |

## Column-ISR / DMA marshaling cost

Same complete windows as the scope trees.

| Preset / ISR | Calls/frame | Min/avg/max µs | CPU share |
|---|--:|--:|--:|
| Shell Flight / isr_wake | 2301.8 | 0.363/1.571/12.041 | 2.89% |
| Shell Flight / isr_pack | 287.7 | 6.006/6.713/9.418 | 1.54% |
| Shell Flight / isr_dma_submit | 287.7 | 0.700/0.936/1.046 | 0.21% |
| Shell Close Flight / isr_wake | 2307.9 | 0.363/1.593/13.241 | 2.93% |
| Shell Close Flight / isr_pack | 288.5 | 6.015/6.871/10.731 | 1.58% |
| Shell Close Flight / isr_dma_submit | 288.5 | 0.708/0.936/1.020 | 0.21% |

Packing marshals LED data on the CPU; submission launches asynchronous DMA. A 600-byte transfer takes approximately 400 µs at 12 MHz. ISR time is already included in render scopes.

ISR CPU shares total 4.64–4.72%, leaving about 59.5 ms foreground time per display window. Peak render needs approximately 1.88× improvement for Shell Flight and 1.94× for Shell Close Flight to fit 62.5 ms.

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

## Global -O3 vs selective -O3

| Preset | Shipping mean ms | O3 mean ms | Mean speedup | FLASH code Δ | RAM1 code Δ |
|---|--:|--:|--:|--:|--:|
| Shell Flight | 121.418 | 115.243 | 1.054× | +13,352 B | +4,736 B |
| Shell Close Flight | 120.839 | 115.208 | 1.049× | +13,352 B | +4,736 B |
