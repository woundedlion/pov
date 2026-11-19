# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\rsjerk-20261118\qualified-color-com3.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `1049347f77b3afa286ab3be6ad28edd0ed1ee884` |
| Method | Show scope = 105 s from profile instance construction. Captured frames 1 through 1668 within scope. Counter window 32. Raw capture 142 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh RiverStones profile 109 32`; analyze only the 105 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:125576, data:156732, headers:8508; free:1740800 / RAM1: variables:314080, code:34856, padding:30680; free:144672 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1101816, data:430372, headers:8932; free:490496 / RAM1: variables:313760, code:191512, padding:5096; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 961–992 root cycles / 600 MHz = 1815445.138 us; measured wall sum 1815444 us, delta 0.63 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 23:32, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `3f3fed2e4c54db98e9823e9a55a7625971f4048d8e0a5c8099cdf13630eaaa42`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Both preset transitions

The 105-second opener covers both 600-frame dwells and the 480-frame forward morph. The reverse morph runs in continuous playback. The table uses both complete 142-second captures. Frame 1 is setup.

| Board | Phase | Frames | Mean (ms) | p99 (ms) | Peak (ms) | Frames at or above 59 ms |
|---|---|---:|---:|---:|---:|---:|
| COM3 | Preset 1 | 598 | 48.985 | 54.059 | 54.758 | 0 |
| COM3 | 1 to 2 | 480 | 54.218 | 57.851 | 58.167 | 0 |
| COM3 | Preset 2 | 599 | 48.328 | 54.106 | 56.041 | 0 |
| COM3 | 2 to 1 | 480 | 54.256 | 58.231 | 58.811 | 0 |
| COM3 | All live frames | 2247 | 51.106 | 57.799 | 58.811 | 0 |
| COM4 | Preset 1 | 598 | 48.987 | 54.030 | 54.847 | 0 |
| COM4 | 1 to 2 | 480 | 54.215 | 57.863 | 58.180 | 0 |
| COM4 | Preset 2 | 599 | 48.327 | 54.039 | 56.087 | 0 |
| COM4 | 2 to 1 | 480 | 54.250 | 58.153 | 58.805 | 0 |
| COM4 | All live frames | 2247 | 51.104 | 57.833 | 58.805 | 0 |

Both boards kept all live frames below 59 ms. The shipping-image check on COM4 ran for 113 seconds and crossed the first handoff with no epoch trap. All runs used PJRC Teensy Loader CLI 2.3.

| Board | Extended capture | SHA-256 |
|---|---|---|
| COM3 | `qualified-color-com3.log` | `3f3fed2e4c54db98e9823e9a55a7625971f4048d8e0a5c8099cdf13630eaaa42` |
| COM4 | `periodic-color-com4.log` | `4bb265476584f1d082b6ab71452c064aa14e2ccc8d411b21c5dcb6e1b76d7abd` |

### Recipe changes and fidelity

Both hex fields follow the same clock using their respective lattice periods. During the source morph, the hue layer evaluates the direct noise and the displacement using their authored strengths. It interpolates the unwrapped hue turns before the hue, lightness, and opacity finish. Daydream samples the hex field analytically.

On Teensy, the two hex fields use precomputed periodic triangles during an authored blend. The grid data occupy 13,068 bytes. Settled presets matched the exact RGB16 comparison. Across 97 sampled frames per direction at 288 by 144 in the quadrant renderer, the periodic sampler differed from the analytic renderer by 12.41 and 14.23 RGB16 RMS codes. The corresponding 99th-percentile absolute differences were 54 and 65 out of 65,535. Independent source tests measured maximum absolute errors of 0.01013 and 0.00366 in signed field units.

The two-preset scroll correction before the hue change peaked at 58.721 ms. With the hue blend and the earlier sector lookup, the peak was 59.963 ms. The earlier 58.703 ms result covered the one-preset, 30-second opener.

## Frame cadence

The playlist assigns **105 s**. Show-frame peak render is **70.027 ms**, spilled **1/1668**. Setup frame 1 is **70.027 ms** and remains separate. Live show peak is **58.167 ms**, spilled **0/1667**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1668. The raw log retains 2248 frame records; the 580 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 58.167 (2), spilled 🟢 0/1667 (0.00%).

README setup frame 1: 70.027 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    43.239 ms  25.944 Mcyc 100.0%
  pov_preserve_half      134.91 us   80.95 kcyc   0.3% x0.96875 139 us/call
  fx_shader_draw         37.551 ms  22.530 Mcyc  86.8%
    fx_prepare_column     1.271 ms  762.76 kcyc   2.9% x144 9 us/call
  fx_prepare_frame        5.133 ms   3.080 Mcyc  11.9% x1 5133 us/call
  fx_advance              64.37 us   38.62 kcyc   0.1% x1 64 us/call
  fx_timeline_step        73.07 us   43.84 kcyc   0.2% x1 73 us/call
  canvas_clear            95.04 us   57.03 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait       1.34 us    0.80 kcyc   0.0% x1 1 us/call
```

Wall min/avg/max = 40.433/43.239/70.027 ms. `fx_shader_draw` averages 37.551 ms/f; total render averages 43.238 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 961–992)

```
frame                    56.733 ms  34.040 Mcyc 100.0%
  pov_preserve_half      136.22 us   81.73 kcyc   0.2% x1 136 us/call
  fx_shader_draw         50.554 ms  30.332 Mcyc  89.1%
    fx_prepare_column     2.204 ms   1.322 Mcyc   3.9% x144 15 us/call
  fx_prepare_frame        5.410 ms   3.246 Mcyc   9.5% x1 5410 us/call
  fx_advance              82.13 us   49.28 kcyc   0.1% x1 82 us/call
  fx_timeline_step       242.80 us  145.68 kcyc   0.4% x1 243 us/call
  canvas_clear            93.50 us   56.10 kcyc   0.2% x1 93 us/call
  canvas_buffer_wait       1.49 us    0.89 kcyc   0.0% x1 1 us/call
```

Wall min/avg/max = 55.812/56.732/58.168 ms. `fx_shader_draw` averages 50.554 ms/f; total render averages 56.731 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 193–224)

```
frame                    52.986 ms  31.792 Mcyc 100.0%
  pov_preserve_half      136.28 us   81.77 kcyc   0.3% x1 136 us/call
  fx_shader_draw         47.025 ms  28.215 Mcyc  88.7%
    fx_prepare_column     1.355 ms  813.15 kcyc   2.6% x144 9 us/call
  fx_prepare_frame        5.395 ms   3.237 Mcyc  10.2% x1 5395 us/call
  fx_advance              65.67 us   39.40 kcyc   0.1% x1 66 us/call
  fx_timeline_step        79.30 us   47.58 kcyc   0.1% x1 79 us/call
  canvas_clear            91.33 us   54.80 kcyc   0.2% x1 91 us/call
  canvas_buffer_wait       1.33 us    0.80 kcyc   0.0% x1 1 us/call
```

Wall min/avg/max = 51.176/52.986/54.759 ms. `fx_shader_draw` averages 47.025 ms/f; total render averages 52.985 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 993–1024)

```
frame                    57.088 ms  34.253 Mcyc 100.0%
  pov_preserve_half      136.48 us   81.89 kcyc   0.2% x1 136 us/call
  fx_shader_draw         50.923 ms  30.554 Mcyc  89.2%
    fx_prepare_column     2.219 ms   1.332 Mcyc   3.9% x144 15 us/call
  fx_prepare_frame        5.398 ms   3.239 Mcyc   9.5% x1 5398 us/call
  fx_advance              81.66 us   49.00 kcyc   0.1% x1 82 us/call
  fx_timeline_step       242.63 us  145.58 kcyc   0.4% x1 243 us/call
  canvas_clear            91.03 us   54.62 kcyc   0.2% x1 91 us/call
  canvas_buffer_wait       1.48 us    0.89 kcyc   0.0% x1 1 us/call
```

Wall min/avg/max = 55.137/57.087/58.105 ms. `fx_shader_draw` averages 50.923 ms/f; total render averages 57.086 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1633–1664)

```
frame                    47.307 ms  28.384 Mcyc 100.0%
  pov_preserve_half      137.63 us   82.58 kcyc   0.3% x1 138 us/call
  fx_shader_draw         41.360 ms  24.816 Mcyc  87.4%
    fx_prepare_column     1.327 ms  796.32 kcyc   2.8% x144 9 us/call
  fx_prepare_frame        5.385 ms   3.231 Mcyc  11.4% x1 5385 us/call
  fx_advance              59.77 us   35.86 kcyc   0.1% x1 60 us/call
  fx_timeline_step        75.40 us   45.24 kcyc   0.2% x1 75 us/call
  canvas_clear            90.91 us   54.55 kcyc   0.2% x1 91 us/call
  canvas_buffer_wait       1.56 us    0.94 kcyc   0.0% x1 2 us/call
```

Wall min/avg/max = 46.209/47.306/50.419 ms. `fx_shader_draw` averages 41.360 ms/f; total render averages 47.305 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

The 105-second opener includes both dwells and the forward morph. The extended capture also covers the reverse morph in continuous playback.

The show-scoped data covers 2 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · RiverStones | 18 | 193–224 | held | — | — | 47.025 | 52.985 | 8 | 54.758 | 0/598 |
| 2 · RiverStones 2 | 18 | 1153–1184 | held | — | — | 45.392 | 51.339 | 16 | 58.167 | 0/1069 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 961–992; CPU shares use its separately snapshotted ISR interval of 2003391 us.

```
isr_wake       1154.0/f 0.471/1.781/28.713 us 3.28%
  isr_pack       144.2/f 6.453/7.570/10.371 us 1.74%
  isr_dma_submit 144.2/f 0.620/0.910/4.060 us 0.20%
```

- Mean pack 7.570 us/call versus submit 0.910 us/call (8.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.28%, corresponding to 60.450 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.692×; Peak live show render: 0.908×; Costliest qualified preset: 0.848×; Transition / non-hold regime: 0.913×; Show end: 0.757×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 89.1% of the frame, 50.554 ms/f.
2. `fx_prepare_frame` — 9.5% of the frame, 5.410 ms/f.
3. `fx_prepare_column` — 3.9% of the frame, 2.204 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 1668 completed by 104.981267 s from profile instance construction. The full-playlist handoff is unmeasured.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Historical single-preset qualification

RiverStones opens the show for 30 seconds. Both desktop Teensys ran the revised effect with PJRC Teensy Loader CLI 2.3. The measured source tree is identical to the integration commit `53056f51c13a1e9c52b54f83b843b1f9e601d5d2`.

| Teensy | Live frames | Mean (ms) | p99 (ms) | Peak (ms) | Frames at or above 59 ms |
|---|---:|---:|---:|---:|---:|
| COM4 | 472 | 52.540 | 58.366 | 58.703 | 0 |
| COM3 | 472 | 52.539 | 58.417 | 58.720 | 0 |

The six native RGB16 captures match the pre-review shipping output for all 480 frames. The RiverStones digest is `a0ed54cd5dacf74c9fcdaeef1796e5ae7c7491ae6ed1f8323c5b6065d1138251`. Both arms now share one bound evaluator per column time. The phase-sampled effects share their endpoint fallback. The stationary simplex cache can serve any admitted color mechanic. Packing its six 23 by 23 cell grids removes 5,640 unused bytes before key bookkeeping.

| Shipping flash (bytes) | Before review | After review |
|---|---:|---:|
| Full show code | 1,092,696 | 1,089,024 |
| Full show data | 429,268 | 429,272 |
| Marginal RiverStones code | 16,000 | 15,800 |
| Marginal RiverStones data | 488 | 496 |

Each marginal value subtracts a build of the same source with only the RiverStones playlist entry excluded.

The general instruction loop peaked at 65.604 ms on the opener. Compiling independent layers from the preset and reusing identical byte-cube corners put both boards below 59 ms. Sharing the palette finish brought marginal code to 15,800 bytes. An alternative that size-optimized stationary cache setup and inlined the palette finish peaked at 58.228 ms across both boards. It used 16,072 marginal code bytes.

Shipping RAM1 code grew from 185,576 to 190,888 bytes within the same six instruction banks. Padding fell from 11,032 to 5,720 bytes. RAM1 variables occupy 313,760 bytes, leaving 13,920 bytes for the stack. RAM2 variables occupy 520,064 bytes, leaving 4,224 bytes free.

## Earlier artwork and flash qualification

The following measurements describe the pre-review shipping build at `b8adc221c`.


RiverStones opens the shipping show for 30 seconds. The revised desktop preset uses a 1/64 burst rate, a translating affine warp on the hexagonal field, and the outer camera with zero spin and 0.1 wander. The measurements below include all live frames within the opening slot on each board. Compared with the previously installed two-preset renderer, shipping FLASH code rose by 568 bytes and FLASH data fell by 928 bytes. RAM1 code fell by 144 bytes and RAM1 variables stayed at 313,760 bytes. In the same-board 35-second Iridescence check, peak render went from 45.154 to 45.326 ms. Both runs had zero display misses.

Both boards used PJRC Teensy Loader CLI 2.3.

| Teensy | Live frames | Mean (ms) | p99 (ms) | Peak (ms) | Frames at or above 59 ms |
|---|---:|---:|---:|---:|---:|
| COM4 | 472 | 52.685 | 58.568 | 58.922 | 0 |
| COM3 | 472 | 52.682 | 58.550 | 58.868 | 0 |

Two additional captures peaked at 58.879 ms on COM4 and 58.865 ms on COM3. All 1,888 live frames across the four captures finished below 59 ms.

The previous shared renderer peaked at 75.015 ms on COM4 and missed the 62.5 ms display deadline on 381 of 433 live frames. Its handwritten predecessor peaked at 79.542 ms with the same revised artwork. The driver primes both halves of the display before it starts the usual quadrant updates. That first draw is recorded separately below.

The ripple evaluator groups waves that share a center, skips disjoint support intervals and preserves the analytic displacement in launch order. The two animated byte cubes share symmetric sampling directions and evaluate each seam once. The three color layers reuse cube coordinates. Only the stationary hue noise uses a spatial approximation. Eligible instances fit a cubic patch in each cube cell. Edits to the noise recipe or generator that invalidate this cache use the direct evaluator. The affine motion, color mechanics and camera remain active.

Across 480 native quadrant frames, the direct-noise control reproduced the reference RGB16 file exactly (SHA-256 `4435faa3cdf266170785fb49847511d65f0dbbfc21b8af9da8c2fe22b00679ad`). The shipped cubic sampler matched 3,610,635 of 4,976,640 pixels exactly. Mean absolute channel error was 0.260 out of 65,535. The 99th percentile of the worst channel per pixel was 6, the maximum was 380, and 12 pixels differed by more than 257 in any channel (one linear 8-bit step). The shipped output SHA-256 was `a0ed54cd5dacf74c9fcdaeef1796e5ae7c7491ae6ed1f8323c5b6065d1138251`.

| Shipping flash (bytes) | Previous renderer | This build |
|---|---:|---:|
| Full show code | 1,172,064 | 1,092,696 |
| Full show data | 430,244 | 429,268 |
| Marginal RiverStones code | 16,664 | 16,000 |
| Marginal RiverStones data | 404 | 488 |

Each marginal figure subtracts a shipping build of the same source with only the RiverStones playlist row excluded. The fitted noise coefficients live in the effect arena. They are generated when needed and require no coefficient array in flash.

### Earlier shared kernels and flash attribution

A frame-interval sampler prepares one column time from the two frame endpoints. Opposing columns use the same sweep time. Each distinct frame layout gets a typed sampler. The first sub-59 ms build contained 19 sampler bodies totaling 29,924 bytes. This build has 19 totaling 14,136 bytes. Counts come from unique linked addresses, including three symbols the default ARM demangler leaves encoded.

The first sub-59 ms build added 21,192 bytes of marginal code, up 4,528 from 16,664. In matched shipping ELFs the hue/lift/mask shader grew by 1,234 bytes, the new ripple evaluator added 1,244, the ripple bounds grew by 584, and the cubic noise routines added 1,230. RiverStones' 1,940-byte column-state sampler was already present before the increase. The original cubic build and this revision produce the same 480-frame RGB16 digest. The pixel-identical direct-noise control peaked at 64.505 ms before this flash pass.

The flash pass reduced RiverStones from 21,192 to 16,000 marginal code bytes and from 500 to 488 marginal data bytes. Its shipping frame sampler fell from 1,940 to 1,164 bytes, the dedicated color shader from 2,780 to 1,718 bytes, and the supported-ripple evaluator from 1,244 to 752 bytes. The palette now records uniform alpha when its table changes, eliminating 70,078 repeated 256-entry scans in the 480-frame native run.

The initial flash-optimized build raised the Iridescence peak to 48.252 ms. Placing its phase-shading and byte-cube routines in instruction RAM recovered the time. The final Iridescence qualification on COM3 peaked at 45.515 ms, compared with 47.242 ms before the RiverStones flash pass. On the full shipping build the placement moved 2,032 additional bytes into instruction RAM, using linker padding. The free RAM1 space remained 13,920 bytes. The RiverStones profile ELF had byte-identical loadable sections before and after this placement. Its earlier 10 kHz program-counter capture put 19.2% of render samples in the hue/lift/mask shader, 18.3% in the supported-ripple evaluator and 8.0% in the shared perceptual color finish. Iridescence, AshCloud, GridWarp, Involution and Millefiori matched their reference RGB16 captures.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
