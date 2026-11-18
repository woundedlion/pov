# RiverStones on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river59-evidence\orientation4.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RiverStones 288×144, single-entry playlist, source `9231cd43a8d33d4d59341dbce3541bce8ed60524` |
| Method | Show scope = 30 s from profile instance construction. Captured frames 1 through 473 within scope. Counter window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh RiverStones profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:115400, data:156812, headers:8364; free:1751040 / RAM1: variables:314080, code:33384, padding:32152; free:144672 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1164408, data:430344, headers:8816; free:428048 / RAM1: variables:313760, code:182632, padding:13976; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 1994964.643 us; measured wall sum 1994967 us, delta 1.18 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 13:42, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `15dea02052f77c619fd29b7a5071773311f4636c98ef3629da42fc6baa5f6ef9`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Live-frame budget and artwork fidelity

RiverStones opens the shipping show for 30 seconds. The revised desktop preset uses a 1/64 burst rate, a translating affine warp on the hexagonal field, and the outer camera with zero spin and 0.1 wander. The measurements below include all live frames within the opening slot on each board. Both boards used PJRC Teensy Loader CLI 2.3.

| Teensy | Live frames | Mean (ms) | p99 (ms) | Peak (ms) | Frames at or above 59 ms |
|---|---:|---:|---:|---:|---:|
| COM4 | 472 | 52.139 | 57.892 | 58.260 | 0 |
| COM3 | 472 | 52.147 | 58.078 | 58.270 | 0 |

The previous shared renderer peaked at 75.015 ms on COM4 and missed the 62.5 ms display deadline on 381 of 433 live frames. Its handwritten predecessor peaked at 79.542 ms with the same revised artwork. The driver primes both halves of the display before it starts the usual quadrant updates. That first draw is recorded separately below.

The ripple evaluator groups waves that share a center, skips disjoint support intervals and preserves the analytic displacement in launch order. The two animated byte cubes share symmetric sampling directions and evaluate each seam once. The three color layers reuse cube coordinates. Only the stationary hue noise uses a spatial approximation. Eligible instances fit a cubic patch in each cube cell. Edits to the noise recipe or generator that invalidate this cache use the direct evaluator. The affine motion, color mechanics and camera remain active.

Across 480 native quadrant frames, the direct-noise control reproduced the reference RGB16 file exactly (SHA-256 `4435faa3cdf266170785fb49847511d65f0dbbfc21b8af9da8c2fe22b00679ad`). The shipped cubic sampler matched 3,610,635 of 4,976,640 pixels exactly. Mean absolute channel error was 0.260 out of 65,535. The 99th percentile of the worst channel per pixel was 6, the maximum was 380, and 12 pixels differed by more than 257 in any channel (one linear 8-bit step). The shipped output SHA-256 was `a0ed54cd5dacf74c9fcdaeef1796e5ae7c7491ae6ed1f8323c5b6065d1138251`.

| Shipping flash (bytes) | Previous renderer | This build |
|---|---:|---:|
| Full show code | 1,172,064 | 1,164,408 |
| Full show data | 430,244 | 430,344 |
| Marginal RiverStones code | 16,664 | 21,192 |
| Marginal RiverStones data | 404 | 500 |

Each marginal figure subtracts a shipping build of the same source with only the RiverStones playlist row excluded. The fitted noise coefficients live in the effect arena. They are generated when needed and require no coefficient array in flash.

The pixel-identical direct-noise control peaked at 64.505 ms. Iridescence, which also uses the shared renderer, peaked at 46.684 ms across 479 live frames on COM3 and missed no display windows.

A separate 10 kHz program-counter capture collected 214,598 render samples on COM4. The combined hue, lift and mask shader accounted for 27.4%, ripple displacement for 24.5%, and sphere-noise lookup table preparation for 8.7%. No samples were missed or dropped. In the shipping ELF, excluding RiverStones removed 2,780 bytes of ITCM for that shader, 1,940 bytes of flash-resident code for frame-interval evaluation, and 1,244 bytes of ITCM for ripple displacement.

## Frame cadence

The playlist assigns **30 s**. Show-frame peak render is **69.253 ms**, spilled **1/473**. Setup frame 1 is **69.253 ms** and remains separate. Live show peak is **58.260 ms**, spilled **0/472**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–473. The raw log retains 520 frame records; the 47 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 58.260, spilled 🟢 0/472 (0.00%).

README setup frame 1: 69.253 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.045 ms  37.227 Mcyc 100.0%
  pov_preserve_half      133.68 us   80.21 kcyc   0.2% x0.96875 138 us/call
  fx_shader_draw         39.063 ms  23.438 Mcyc  63.0%
    fx_prepare_column     1.709 ms   1.025 Mcyc   2.8% x144 12 us/call
  fx_prepare_frame        4.669 ms   2.801 Mcyc   7.5% x1 4669 us/call
  fx_advance              53.35 us   32.01 kcyc   0.1% x1 53 us/call
  fx_timeline_step        81.12 us   48.67 kcyc   0.1% x1 81 us/call
  canvas_clear            91.89 us   55.13 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     17.764 ms  10.658 Mcyc  28.6% x1 17764 us/call
```

Wall min/avg/max = 41.191/62.044/69.253 ms. `fx_shader_draw` averages 39.063 ms/f; total render averages 44.281 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.343 ms  37.406 Mcyc 100.0%
  pov_preserve_half      136.88 us   82.13 kcyc   0.2% x1 137 us/call
  fx_shader_draw         51.220 ms  30.732 Mcyc  82.2%
    fx_prepare_column     2.077 ms   1.246 Mcyc   3.3% x144 14 us/call
  fx_prepare_frame        4.674 ms   2.804 Mcyc   7.5% x1 4674 us/call
  fx_advance              52.34 us   31.40 kcyc   0.1% x1 52 us/call
  fx_timeline_step        77.40 us   46.44 kcyc   0.1% x1 77 us/call
  canvas_clear            87.79 us   52.68 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait      5.913 ms   3.548 Mcyc   9.5% x1 5913 us/call
```

Wall min/avg/max = 59.981/62.342/65.625 ms. `fx_shader_draw` averages 51.220 ms/f; total render averages 56.430 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 417–448)

```
frame                    62.250 ms  37.350 Mcyc 100.0%
  pov_preserve_half      136.63 us   81.98 kcyc   0.2% x1 137 us/call
  fx_shader_draw         49.936 ms  29.962 Mcyc  80.2%
    fx_prepare_column     2.029 ms   1.217 Mcyc   3.3% x144 14 us/call
  fx_prepare_frame        4.669 ms   2.802 Mcyc   7.5% x1 4669 us/call
  fx_advance              53.54 us   32.12 kcyc   0.1% x1 54 us/call
  fx_timeline_step        77.39 us   46.44 kcyc   0.1% x1 77 us/call
  canvas_clear            87.86 us   52.71 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait      7.106 ms   4.264 Mcyc  11.4% x1 7106 us/call
```

Wall min/avg/max = 60.454/62.249/64.580 ms. `fx_shader_draw` averages 49.936 ms/f; total render averages 55.144 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 1999475 us.

```
isr_wake       1151.8/f 0.530/1.748/23.700 us 3.22%
  isr_pack       144.0/f 6.556/7.410/10.245 us 1.70%
  isr_dma_submit 144.0/f 0.598/0.915/1.650 us 0.21%
```

- Mean pack 7.410 us/call versus submit 0.915 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.22%, corresponding to 60.487 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.708×; Peak live show render: 0.903×; Show end: 0.882×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 82.2% of the frame, 51.220 ms/f.
2. `fx_prepare_frame` — 7.5% of the frame, 4.674 ms/f.
3. `fx_prepare_column` — 3.3% of the frame, 2.077 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- The show cutoff uses the first elapsed counter window, subsequent ISR snapshot intervals, and per-frame wall times in the boundary window. Frame 473 completed by 29.987835 s from profile instance construction. The full-playlist handoff is unmeasured.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RiverStones`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
