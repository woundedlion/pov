# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river59-evidence\review-finish4.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `62520ec25e068949eeb144e34cfd023ef05f4daf` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1 through 473 within scope. Counter window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:111224, data:156424, headers:8832; free:1755136 / RAM1: variables:314080, code:33816, padding:31720; free:144672 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1089024, data:429272, headers:8480; free:504840 / RAM1: variables:313760, code:190888, padding:5720; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 1996665.317 us; measured wall sum 1996666 us, delta 0.34 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 18:11, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `6d026853f0f1a04adcd65496346a754da0b1ed46d515e1e542ab6254a1814157`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Review qualification

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


RiverStones opens the shipping show for 30 seconds. The revised desktop preset uses a 1/64 burst rate, a translating affine warp on the hexagonal field, and the outer camera with zero spin and 0.1 wander. The measurements below include all live frames within the opening slot on each board. Both boards used PJRC Teensy Loader CLI 2.3.

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

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **71.753 ms**, spilled **1/473**. Setup frame 1 is **71.753 ms** and remains separate. Live show peak is **58.703 ms**, spilled **0/472**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–473. The raw log retains 520 frame records; the 47 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 58.703, spilled 🟢 0/472 (0.00%).

README setup frame 1: 71.753 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.187 ms  37.312 Mcyc 100.0%
  pov_preserve_half      134.83 us   80.90 kcyc   0.2% x0.96875 139 us/call
  fx_shader_draw         39.236 ms  23.542 Mcyc  63.1%
    fx_prepare_column     1.253 ms  751.83 kcyc   2.0% x144 9 us/call
  fx_prepare_frame        4.892 ms   2.935 Mcyc   7.9% x1 4892 us/call
  fx_advance              52.77 us   31.66 kcyc   0.1% x1 53 us/call
  fx_timeline_step        70.00 us   42.00 kcyc   0.1% x1 70 us/call
  canvas_clear            91.05 us   54.63 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     17.538 ms  10.523 Mcyc  28.2% x1 17538 us/call
```

Wall min/avg/max = 41.411/62.186/71.755 ms. `fx_shader_draw` averages 39.236 ms/f; total render averages 44.648 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.396 ms  37.437 Mcyc 100.0%
  pov_preserve_half      137.28 us   82.37 kcyc   0.2% x1 137 us/call
  fx_shader_draw         51.506 ms  30.903 Mcyc  82.5%
    fx_prepare_column     1.566 ms  939.58 kcyc   2.5% x144 11 us/call
  fx_prepare_frame        4.895 ms   2.937 Mcyc   7.8% x1 4895 us/call
  fx_advance              54.07 us   32.44 kcyc   0.1% x1 54 us/call
  fx_timeline_step        68.39 us   41.03 kcyc   0.1% x1 68 us/call
  canvas_clear            87.10 us   52.26 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait      5.482 ms   3.289 Mcyc   8.8% x1 5482 us/call
```

Wall min/avg/max = 60.045/62.395/65.881 ms. `fx_shader_draw` averages 51.506 ms/f; total render averages 56.913 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 417–448)

```
frame                    62.298 ms  37.379 Mcyc 100.0%
  pov_preserve_half      137.44 us   82.47 kcyc   0.2% x1 137 us/call
  fx_shader_draw         50.215 ms  30.129 Mcyc  80.6%
    fx_prepare_column     1.541 ms  924.83 kcyc   2.5% x144 11 us/call
  fx_prepare_frame        4.893 ms   2.936 Mcyc   7.9% x1 4893 us/call
  fx_advance              52.99 us   31.79 kcyc   0.1% x1 53 us/call
  fx_timeline_step        68.78 us   41.27 kcyc   0.1% x1 69 us/call
  canvas_clear            86.98 us   52.19 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait      6.678 ms   4.007 Mcyc  10.7% x1 6678 us/call
```

Wall min/avg/max = 60.545/62.298/64.539 ms. `fx_shader_draw` averages 50.215 ms/f; total render averages 55.620 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 1999471 us.

```
isr_wake       1151.8/f 0.500/1.743/24.831 us 3.21%
  isr_pack       144.0/f 6.526/7.400/10.591 us 1.70%
  isr_dma_submit 144.0/f 0.593/0.901/1.256 us 0.20%
```

- Mean pack 7.400 us/call versus submit 0.901 us/call (8.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.21%, corresponding to 60.494 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.714×; Peak live show render: 0.911×; Show end: 0.890×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 82.5% of the frame, 51.506 ms/f.
2. `fx_prepare_frame` — 7.8% of the frame, 4.895 ms/f.
3. `fx_prepare_column` — 2.5% of the frame, 1.566 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 473 completed by 29.987932 s from profile instance construction. The full-playlist handoff is unmeasured.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
