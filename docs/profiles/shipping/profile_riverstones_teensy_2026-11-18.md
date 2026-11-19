# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

The ranked results cover the 30-second show slot. The extended capture also covers both preset transitions.
Raw capture: `C:\work\temp\rsjerk-20261118\production-com3.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `a6ae35c5e45070518210a6ebb725309a7827a7cd` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1 through 473 within scope. Counter window 32. Raw capture 142 s; later frames appear in the extended cycle measurements. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:124624, data:156732, headers:8436; free:1741824 / RAM1: variables:314080, code:34392, padding:31144; free:144672 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1099248, data:430372, headers:8404; free:493592 / RAM1: variables:313760, code:190408, padding:6200; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 1996515.615 us; measured wall sum 1996522 us, delta 3.20 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 22:37, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `1b117734ceb784b1cbb4a69316edb0c621f5510098eacba2c98531f7c88d0fa4`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Two-preset motion and timing

The 30-second opener contains only preset 1. Each 142-second capture continued into both 480-frame transitions. The peak across all captured live frames on both boards was **58.721 ms**. Frame 1 is setup. The show ranking uses only its 30-second slot.

| Board | Phase | Live frames | Mean (ms) | p99 (ms) | Peak (ms) | At or above 59 ms | Above 62.5 ms |
|---|---|---:|---:|---:|---:|---:|---:|
| COM3 | Preset 1 | 598 | 49.189 | 54.384 | 55.151 | 0 | 0 |
| COM3 | 1 to 2 | 480 | 52.035 | 56.124 | 56.652 | 0 | 0 |
| COM3 | Preset 2 | 599 | 48.304 | 54.056 | 55.923 | 0 | 0 |
| COM3 | 2 to 1 | 480 | 54.119 | 58.310 | 58.721 | 0 | 0 |
| COM3 | Entire capture after setup | 2247 | 50.668 | 57.413 | 58.721 | 0 | 0 |
| COM4 | Preset 1 | 598 | 49.188 | 54.241 | 55.287 | 0 | 0 |
| COM4 | 1 to 2 | 480 | 52.039 | 56.117 | 56.636 | 0 | 0 |
| COM4 | Preset 2 | 599 | 48.307 | 54.075 | 56.021 | 0 | 0 |
| COM4 | 2 to 1 | 480 | 54.120 | 58.267 | 58.673 | 0 | 0 |
| COM4 | Entire capture after setup | 2247 | 50.669 | 57.377 | 58.673 | 0 | 0 |

The in-show capture on COM3 includes frames 2 through 473. Setup took 70.523 ms on COM3 and 70.426 ms on COM4. The first destination marker precedes frame 600 and the return marker precedes frame 1679. The table assigns the next 480 frames to each blend, then resumes the destination dwell.

### Before the correction

The previous source-value morph used one affine offset for two different hex-cell sizes. When the scroll clock wrapped every 50 frames, one partially visible field moved by a fraction of a cell. Both fields now use the same column clock with offsets expressed in their own lattice periods. Their scalar values are blended with the authored easing before color modulation.

| Build | Board | Preset 1 peak | Forward peak | Preset 2 peak | Reverse peak | Transition frames above 62.5 ms |
|---|---|---:|---:|---:|---:|---:|
| Installed before the fix | COM4 | 59.090 ms | 63.837 ms | 60.157 ms | 65.445 ms | 94 |
| Corrected motion, analytic Teensy kernel | COM3 | 62.509 ms | 66.996 ms | 62.687 ms | 68.448 ms | 564 |
| Shipping correction | COM3 | 55.151 ms | 56.652 ms | 55.923 ms | 58.721 ms | 0 |
| Shipping correction | COM4 | 55.287 ms | 56.636 ms | 56.021 ms | 58.673 ms | 0 |

Across a 2,162-frame Daydream render at 96 by 20, the largest RGB16 step at a wrap relative to the median step in its eight neighboring frames fell from 2.52 to 0.98 forward and from 3.81 to 1.18 in reverse. Both settled presets matched their previous RGB16 frames. The native coordinate test also covered all 19 clock wraps in both directions.

### Device rendering and memory

Daydream uses the analytic hex field. Teensy precomputes two 61 by 61 tables at startup and uses triangular interpolation only while both authored profiles are visible. Both dwell endpoints use the analytic sampler. A changed radius, shape or softness also selects the analytic sampler. Scale is applied at sample time. The tables occupy 14,884 bytes of grid data in the existing arena. The device effect uses 269,186 of 272,384 persistent bytes, including its stationary noise atlas.

In the native comparison of the device sampling path, the settled images matched. Across 288 by 144 quadrant renders sampled every fifth frame plus boundaries, forward and reverse RGB16 channel RMS errors were 4.97 and 4.62 on a 65,535 scale. The 99th-percentile absolute differences were 19 and 20 counts. The maxima were 1,442 and 582 counts.

Compared with the previously installed two-preset renderer, shipping FLASH code rose by 568 bytes and FLASH data fell by 928 bytes. RAM1 code fell by 144 bytes and RAM1 variables stayed at 313,760 bytes. In the same-board 35-second Iridescence check, peak render went from 45.154 to 45.326 ms. Both runs had zero display misses.

Both boards used PJRC Teensy Loader CLI 2.3. The COM4 capture used the same profile ELF as COM3. The original logs and measured firmware hashes are retained under `C:/work/temp/rsjerk-20261118`.

| Board | Raw log | SHA-256 |
|---|---|---|
| COM3 | `production-com3.log` | `1b117734ceb784b1cbb4a69316edb0c621f5510098eacba2c98531f7c88d0fa4` |
| COM4 | `production-com4.log` | `ca76395bd714ac1466d04d81ff5cbec50b5f81906f6b48b29c347aeb5749ff54` |

The earlier 58.703 ms COM4 result covered 472 live frames of the one-preset opener. The earlier 58.512 ms two-preset result used parameter interpolation. This capture exercises the source-value blend.

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **70.523 ms**, spilled **1/473**. Setup frame 1 is **70.523 ms** and remains separate. Live show peak is **55.151 ms**, spilled **0/472**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–473. The raw log retains 2248 frame records; the 1775 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 55.151 (1), spilled 🟢 0/472 (0.00%).

README setup frame 1: 70.523 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.097 ms  37.258 Mcyc 100.0%
  pov_preserve_half      136.13 us   81.68 kcyc   0.2% x0.96875 141 us/call
  fx_shader_draw         37.919 ms  22.751 Mcyc  61.1%
    fx_prepare_column     1.452 ms  870.97 kcyc   2.3% x144 10 us/call
  fx_prepare_frame        5.085 ms   3.051 Mcyc   8.2% x1 5085 us/call
  fx_advance              59.43 us   35.66 kcyc   0.1% x1 59 us/call
  fx_timeline_step        79.58 us   47.75 kcyc   0.1% x1 80 us/call
  canvas_clear            92.32 us   55.39 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     18.525 ms  11.115 Mcyc  29.8% x1 18525 us/call
```

Wall min/avg/max = 40.880/62.097/70.528 ms. `fx_shader_draw` averages 37.919 ms/f; total render averages 43.573 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.391 ms  37.435 Mcyc 100.0%
  pov_preserve_half      135.50 us   81.30 kcyc   0.2% x1 135 us/call
  fx_shader_draw         47.240 ms  28.344 Mcyc  75.7%
    fx_prepare_column     1.510 ms  905.95 kcyc   2.4% x144 10 us/call
  fx_prepare_frame        5.359 ms   3.215 Mcyc   8.6% x1 5359 us/call
  fx_advance              56.40 us   33.84 kcyc   0.1% x1 56 us/call
  fx_timeline_step       100.61 us   60.36 kcyc   0.2% x1 101 us/call
  canvas_clear            88.80 us   53.28 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait      9.207 ms   5.524 Mcyc  14.8% x1 9207 us/call
```

Wall min/avg/max = 59.431/62.391/65.252 ms. `fx_shader_draw` averages 47.240 ms/f; total render averages 53.185 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 417–448)

```
frame                    62.317 ms  37.390 Mcyc 100.0%
  pov_preserve_half      135.96 us   81.58 kcyc   0.2% x1 136 us/call
  fx_shader_draw         44.963 ms  26.978 Mcyc  72.2%
    fx_prepare_column     1.478 ms  886.53 kcyc   2.4% x144 10 us/call
  fx_prepare_frame        5.359 ms   3.215 Mcyc   8.6% x1 5359 us/call
  fx_advance              60.93 us   36.56 kcyc   0.1% x1 61 us/call
  fx_timeline_step        92.16 us   55.30 kcyc   0.1% x1 92 us/call
  canvas_clear            89.20 us   53.52 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     11.409 ms   6.845 Mcyc  18.3% x1 11409 us/call
```

Wall min/avg/max = 60.448/62.316/63.873 ms. `fx_shader_draw` averages 44.963 ms/f; total render averages 50.908 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Only the first preset plays during the 30-second opener. The extended capture below follows both preset holds and both transitions.

The show-scoped data covers 1 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · RiverStones | 14 | 193–224 | held | — | — | 47.240 | 53.185 | 8 | 55.151 | 0/472 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 2000050 us.

```
isr_wake       1152.1/f 0.471/1.740/22.133 us 3.20%
  isr_pack       144.0/f 6.526/7.400/9.995 us 1.70%
  isr_dma_submit 144.0/f 0.583/0.911/8.538 us 0.21%
```

- Mean pack 7.400 us/call versus submit 0.911 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.20%, corresponding to 60.500 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.697×; Peak live show render: 0.851×; Show end: 0.815×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 75.7% of the frame, 47.240 ms/f.
2. `fx_prepare_frame` — 8.6% of the frame, 5.359 ms/f.
3. `fx_prepare_column` — 2.4% of the frame, 1.510 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 473 completed by 29.984595 s from profile instance construction. The full-playlist handoff is unmeasured.
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
