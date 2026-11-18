# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river59-evidence\size-alpha4.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `6f858a9d3f173cc1a5b6b561847443dbd75156ef` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1 through 473 within scope. Counter window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:110224, data:156420, headers:8812; free:1756160 / RAM1: variables:314080, code:31560, padding:1208; free:177440 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1092712, data:429268, headers:8884; free:500752 / RAM1: variables:313760, code:183544, padding:13064; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 1996520.527 us; measured wall sum 1996520 us, delta 0.26 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 14:55, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `9f12808e37c4b9370e1fc7dcbff129097a7511ff06d7d05677e74da10f319872`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Live-frame budget and artwork fidelity

RiverStones opens the shipping show for 30 seconds. The revised desktop preset uses a 1/64 burst rate, a translating affine warp on the hexagonal field, and the outer camera with zero spin and 0.1 wander. The measurements below include all live frames within the opening slot on each board. Both boards used PJRC Teensy Loader CLI 2.3.

| Teensy | Live frames | Mean (ms) | p99 (ms) | Peak (ms) | Frames at or above 59 ms |
|---|---:|---:|---:|---:|---:|
| COM4 | 472 | 52.681 | 58.513 | 58.932 | 0 |
| COM3 | 472 | 52.681 | 58.540 | 58.925 | 0 |

Repeat 30-second captures peaked at 58.951 ms on COM4 and 58.909 ms on COM3. All 1,888 live frames in the four captures finished below 59 ms.

The previous shared renderer peaked at 75.015 ms on COM4 and missed the 62.5 ms display deadline on 381 of 433 live frames. Its handwritten predecessor peaked at 79.542 ms with the same revised artwork. The driver primes both halves of the display before it starts the usual quadrant updates. That first draw is recorded separately below.

The ripple evaluator groups waves that share a center, skips disjoint support intervals and preserves the analytic displacement in launch order. The two animated byte cubes share symmetric sampling directions and evaluate each seam once. The three color layers reuse cube coordinates. Only the stationary hue noise uses a spatial approximation. Eligible instances fit a cubic patch in each cube cell. Edits to the noise recipe or generator that invalidate this cache use the direct evaluator. The affine motion, color mechanics and camera remain active.

Across 480 native quadrant frames, the direct-noise control reproduced the reference RGB16 file exactly (SHA-256 `4435faa3cdf266170785fb49847511d65f0dbbfc21b8af9da8c2fe22b00679ad`). The shipped cubic sampler matched 3,610,635 of 4,976,640 pixels exactly. Mean absolute channel error was 0.260 out of 65,535. The 99th percentile of the worst channel per pixel was 6, the maximum was 380, and 12 pixels differed by more than 257 in any channel (one linear 8-bit step). The shipped output SHA-256 was `a0ed54cd5dacf74c9fcdaeef1796e5ae7c7491ae6ed1f8323c5b6065d1138251`.

| Shipping flash (bytes) | Previous renderer | This build |
|---|---:|---:|
| Full show code | 1,172,064 | 1,092,712 |
| Full show data | 430,244 | 429,268 |
| Marginal RiverStones code | 16,664 | 16,000 |
| Marginal RiverStones data | 404 | 488 |

Each marginal figure subtracts a shipping build of the same source with only the RiverStones playlist row excluded. The fitted noise coefficients live in the effect arena. They are generated when needed and require no coefficient array in flash.

## Shared kernels and flash attribution

A frame-interval sampler prepares one column time from the two frame endpoints. Opposing columns use the same sweep time. Each distinct frame layout gets a typed sampler. The first sub-59 ms build contained 19 sampler bodies totaling 29,924 bytes. This build has 19 totaling 14,136 bytes. Counts come from unique linked addresses, including three symbols the default ARM demangler leaves encoded.

The first sub-59 ms build added 21,192 bytes of marginal code, up 4,528 from 16,664. In matched shipping ELFs the hue/lift/mask shader grew by 1,234 bytes, the new ripple evaluator added 1,244, the ripple bounds grew by 584, and the cubic noise routines added 1,230. RiverStones' 1,940-byte column-state sampler was already present before the increase. The original cubic build and this revision produce the same 480-frame RGB16 digest. The pixel-identical direct-noise control peaked at 64.505 ms before this flash pass.

The flash pass reduced the marginal code from 21,192 to 16,000 bytes and marginal data from 500 to 488 bytes. The complete show saved 71,696 code bytes and 1,076 data bytes relative to that first sub-59 ms build. Its RAM1 code increased 912 bytes with the same 13,920 bytes free for local variables. In the shipping ELF, the RiverStones frame sampler fell from 1,940 to 1,164 bytes, the dedicated color shader from 2,780 to 1,718 bytes, and the supported-ripple evaluator from 1,244 to 752 bytes. The stationary cubic coefficients are generated in RAM.

The final 10 kHz program-counter capture on COM4 recorded 217,322 render samples with zero missed or dropped. The hue/lift/mask shader accounted for 19.2%, the supported-ripple evaluator for 18.3%, and the shared perceptual color finish for 8.0%. The palette now records uniform alpha when its table changes. This eliminates 70,078 repeated 256-entry uniformity scans in the 480-frame native run.

Iridescence peaked at 48.252 ms over the first 479 live frames on COM3 with 0 missed display windows. The previous shipping peak on COM3 was 47.242 ms. Iridescence, AshCloud, GridWarp, Involution and Millefiori all matched their reference RGB16 captures.

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **70.477 ms**, spilled **1/473**. Setup frame 1 is **70.477 ms** and remains separate. Live show peak is **58.932 ms**, spilled **0/472**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–473. The raw log retains 520 frame records; the 47 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 58.932, spilled 🟢 0/472 (0.00%).

README setup frame 1: 70.477 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.141 ms  37.285 Mcyc 100.0%
  pov_preserve_half      133.97 us   80.38 kcyc   0.2% x0.96875 138 us/call
  fx_shader_draw         39.539 ms  23.723 Mcyc  63.6%
    fx_prepare_column     1.179 ms  707.63 kcyc   1.9% x144 8 us/call
  fx_prepare_frame        4.623 ms   2.774 Mcyc   7.4% x1 4623 us/call
  fx_advance              60.30 us   36.18 kcyc   0.1% x1 60 us/call
  fx_timeline_step        63.25 us   37.95 kcyc   0.1% x1 63 us/call
  canvas_clear            93.69 us   56.21 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     17.448 ms  10.469 Mcyc  28.1% x1 17448 us/call
```

Wall min/avg/max = 41.465/62.141/70.481 ms. `fx_shader_draw` averages 39.539 ms/f; total render averages 44.693 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.391 ms  37.435 Mcyc 100.0%
  pov_preserve_half      135.00 us   81.00 kcyc   0.2% x1 135 us/call
  fx_shader_draw         51.943 ms  31.166 Mcyc  83.3%
    fx_prepare_column     1.485 ms  891.30 kcyc   2.4% x144 10 us/call
  fx_prepare_frame        4.627 ms   2.776 Mcyc   7.4% x1 4627 us/call
  fx_advance              61.04 us   36.62 kcyc   0.1% x1 61 us/call
  fx_timeline_step        62.55 us   37.53 kcyc   0.1% x1 63 us/call
  canvas_clear            89.60 us   53.76 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait      5.301 ms   3.180 Mcyc   8.5% x1 5301 us/call
```

Wall min/avg/max = 59.986/62.391/65.952 ms. `fx_shader_draw` averages 51.943 ms/f; total render averages 57.091 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 417–448)

```
frame                    62.296 ms  37.378 Mcyc 100.0%
  pov_preserve_half      136.64 us   81.99 kcyc   0.2% x1 137 us/call
  fx_shader_draw         50.573 ms  30.344 Mcyc  81.2%
    fx_prepare_column     1.453 ms  872.02 kcyc   2.3% x144 10 us/call
  fx_prepare_frame        4.629 ms   2.777 Mcyc   7.4% x1 4629 us/call
  fx_advance              62.33 us   37.40 kcyc   0.1% x1 62 us/call
  fx_timeline_step        65.54 us   39.32 kcyc   0.1% x1 66 us/call
  canvas_clear            89.42 us   53.65 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait      6.567 ms   3.940 Mcyc  10.5% x1 6567 us/call
```

Wall min/avg/max = 60.479/62.296/64.681 ms. `fx_shader_draw` averages 50.573 ms/f; total render averages 55.730 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 1999595 us.

```
isr_wake       1151.8/f 0.571/1.756/21.760 us 3.23%
  isr_pack       144.0/f 6.526/7.411/10.060 us 1.70%
  isr_dma_submit 144.0/f 0.641/0.915/1.265 us 0.21%
```

- Mean pack 7.411 us/call versus submit 0.915 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.23%, corresponding to 60.481 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.715×; Peak live show render: 0.913×; Show end: 0.892×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 83.3% of the frame, 51.943 ms/f.
2. `fx_prepare_frame` — 7.4% of the frame, 4.627 ms/f.
3. `fx_prepare_column` — 2.4% of the frame, 1.485 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 473 completed by 29.988012 s from profile instance construction. The full-playlist handoff is unmeasured.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
