# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river59-evidence\recover-river4.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `4bf24181947f354011ca64b75c00d722bc4c40d2` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1 through 473 within scope. Counter window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:110224, data:156420, headers:8812; free:1756160 / RAM1: variables:314080, code:31560, padding:1208; free:177440 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1092696, data:429268, headers:8900; free:500752 / RAM1: variables:313760, code:185576, padding:11032; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 1996581.273 us; measured wall sum 1996584 us, delta 1.37 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 15:16, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `da2e00924eac1ee38f373daec0f8c01c70c7ce4db505e88e8b2b5375388371da`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Live-frame budget and artwork fidelity

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

## Shared kernels and flash attribution

A frame-interval sampler prepares one column time from the two frame endpoints. Opposing columns use the same sweep time. Each distinct frame layout gets a typed sampler. The first sub-59 ms build contained 19 sampler bodies totaling 29,924 bytes. This build has 19 totaling 14,136 bytes. Counts come from unique linked addresses, including three symbols the default ARM demangler leaves encoded.

The first sub-59 ms build added 21,192 bytes of marginal code, up 4,528 from 16,664. In matched shipping ELFs the hue/lift/mask shader grew by 1,234 bytes, the new ripple evaluator added 1,244, the ripple bounds grew by 584, and the cubic noise routines added 1,230. RiverStones' 1,940-byte column-state sampler was already present before the increase. The original cubic build and this revision produce the same 480-frame RGB16 digest. The pixel-identical direct-noise control peaked at 64.505 ms before this flash pass.

The flash pass reduced RiverStones from 21,192 to 16,000 marginal code bytes and from 500 to 488 marginal data bytes. Its shipping frame sampler fell from 1,940 to 1,164 bytes, the dedicated color shader from 2,780 to 1,718 bytes, and the supported-ripple evaluator from 1,244 to 752 bytes. The palette now records uniform alpha when its table changes, eliminating 70,078 repeated 256-entry scans in the 480-frame native run.

The initial flash-optimized build raised the Iridescence peak to 48.252 ms. Placing its phase-shading and byte-cube routines in instruction RAM recovered the time. The final Iridescence qualification on COM3 peaked at 45.515 ms, compared with 47.242 ms before the RiverStones flash pass. On the full shipping build the placement moved 2,032 additional bytes into instruction RAM, using linker padding. The free RAM1 space remained 13,920 bytes. The RiverStones profile ELF had byte-identical loadable sections before and after this placement. Its earlier 10 kHz program-counter capture put 19.2% of render samples in the hue/lift/mask shader, 18.3% in the supported-ripple evaluator and 8.0% in the shared perceptual color finish. Iridescence, AshCloud, GridWarp, Involution and Millefiori matched their reference RGB16 captures.

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **70.453 ms**, spilled **1/473**. Setup frame 1 is **70.453 ms** and remains separate. Live show peak is **58.922 ms**, spilled **0/472**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–473. The raw log retains 520 frame records; the 47 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 58.922, spilled 🟢 0/472 (0.00%).

README setup frame 1: 70.453 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.135 ms  37.281 Mcyc 100.0%
  pov_preserve_half      134.77 us   80.86 kcyc   0.2% x0.96875 139 us/call
  fx_shader_draw         39.538 ms  23.723 Mcyc  63.6%
    fx_prepare_column     1.184 ms  710.14 kcyc   1.9% x144 8 us/call
  fx_prepare_frame        4.625 ms   2.775 Mcyc   7.4% x1 4625 us/call
  fx_advance              63.04 us   37.82 kcyc   0.1% x1 63 us/call
  fx_timeline_step        66.48 us   39.89 kcyc   0.1% x1 66 us/call
  canvas_clear            93.22 us   55.93 kcyc   0.2% x1 93 us/call
  canvas_buffer_wait     17.436 ms  10.462 Mcyc  28.1% x1 17436 us/call
```

Wall min/avg/max = 41.540/62.135/70.457 ms. `fx_shader_draw` averages 39.538 ms/f; total render averages 44.699 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.393 ms  37.436 Mcyc 100.0%
  pov_preserve_half      135.28 us   81.17 kcyc   0.2% x1 135 us/call
  fx_shader_draw         51.941 ms  31.165 Mcyc  83.2%
    fx_prepare_column     1.483 ms  889.85 kcyc   2.4% x144 10 us/call
  fx_prepare_frame        4.627 ms   2.776 Mcyc   7.4% x1 4627 us/call
  fx_advance              62.04 us   37.23 kcyc   0.1% x1 62 us/call
  fx_timeline_step        65.44 us   39.27 kcyc   0.1% x1 65 us/call
  canvas_clear            89.53 us   53.72 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait      5.300 ms   3.180 Mcyc   8.5% x1 5300 us/call
```

Wall min/avg/max = 59.980/62.393/65.931 ms. `fx_shader_draw` averages 51.941 ms/f; total render averages 57.093 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 417–448)

```
frame                    62.296 ms  37.378 Mcyc 100.0%
  pov_preserve_half      136.87 us   82.12 kcyc   0.2% x1 137 us/call
  fx_shader_draw         50.566 ms  30.339 Mcyc  81.2%
    fx_prepare_column     1.449 ms  869.51 kcyc   2.3% x144 10 us/call
  fx_prepare_frame        4.630 ms   2.778 Mcyc   7.4% x1 4630 us/call
  fx_advance              61.60 us   36.96 kcyc   0.1% x1 62 us/call
  fx_timeline_step        63.76 us   38.25 kcyc   0.1% x1 64 us/call
  canvas_clear            89.06 us   53.44 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait      6.570 ms   3.942 Mcyc  10.5% x1 6570 us/call
```

Wall min/avg/max = 60.511/62.296/64.595 ms. `fx_shader_draw` averages 50.566 ms/f; total render averages 55.726 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 1999593 us.

```
isr_wake       1151.8/f 0.571/1.756/17.106 us 3.23%
  isr_pack       144.0/f 6.526/7.415/10.548 us 1.70%
  isr_dma_submit 144.0/f 0.608/0.915/0.981 us 0.21%
```

- Mean pack 7.415 us/call versus submit 0.915 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.23%, corresponding to 60.481 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.715×; Peak live show render: 0.913×; Show end: 0.892×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 83.2% of the frame, 51.941 ms/f.
2. `fx_prepare_frame` — 7.4% of the frame, 4.627 ms/f.
3. `fx_prepare_column` — 2.4% of the frame, 1.483 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 473 completed by 29.987990 s from profile instance construction. The full-playlist handoff is unmeasured.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
