# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\gnomonicstars_ship.log`. Replaces `profile_gnomonicstars_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared scan-shapes, spherical-math and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh GnomonicStars profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:75920, data:148956, headers:8596; free:1798144 / RAM1: variables:314048, code:20872, padding:11896; free:177472 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 2000150.603 us; measured wall sum 2000146 us, delta 2.30 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:40, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `71f6729ce38785db0270f499f66485cca766e53f1097e0066012836cf4c9d4a2`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **26.078 ms**, spilled **0/480**. Setup frame 1 is **15.179 ms** and remains separate. Live show peak is **26.078 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1097 frame records; the 617 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 26.078, spilled 🟢 0/479 (0.00%).

README setup frame 1: 15.179 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    59.361 ms  35.617 Mcyc 100.0%
  pov_preserve_half      139.37 us   83.62 kcyc   0.2% x0.96875 144 us/call
  gn_draw_stars          10.086 ms   6.052 Mcyc  17.0%
    gn_star_scan          8.468 ms   5.081 Mcyc  14.3% x252.812 33 us/call
  gn_spiral_build         17.49 us   10.49 kcyc   0.0% x0.03125 560 us/call
  gn_timeline_step        60.34 us   36.21 kcyc   0.1%
    animation_mobius_step   16.21 us    9.72 kcyc   0.0% x1 16 us/call
  canvas_clear            95.12 us   57.07 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     48.931 ms  29.358 Mcyc  82.4% x1 48931 us/call
```

Wall min/avg/max = 14.758/59.360/72.393 ms. `gn_draw_stars` averages 10.086 ms/f; total render averages 10.431 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.505 ms  37.503 Mcyc 100.0%
  pov_preserve_half      139.88 us   83.93 kcyc   0.2% x1 140 us/call
  gn_draw_stars          18.116 ms  10.870 Mcyc  29.0%
    gn_star_scan         16.602 ms   9.961 Mcyc  26.6% x238.406 70 us/call
  gn_timeline_step        54.40 us   32.64 kcyc   0.1%
    animation_mobius_step   14.28 us    8.57 kcyc   0.0% x1 14 us/call
  canvas_clear            90.91 us   54.55 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     44.073 ms  26.444 Mcyc  70.5% x1 44073 us/call
```

Wall min/avg/max = 48.207/62.504/75.771 ms. `gn_draw_stars` averages 18.116 ms/f; total render averages 18.432 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.457 ms  37.474 Mcyc 100.0%
  pov_preserve_half      145.60 us   87.36 kcyc   0.2% x1 146 us/call
  gn_draw_stars           9.059 ms   5.435 Mcyc  14.5%
    gn_star_scan          7.450 ms   4.470 Mcyc  11.9% x244.406 30 us/call
  gn_timeline_step        57.53 us   34.52 kcyc   0.1%
    animation_mobius_step   17.13 us   10.28 kcyc   0.0% x1 17 us/call
  canvas_clear            91.29 us   54.78 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     53.076 ms  31.845 Mcyc  85.0% x1 53076 us/call
```

Wall min/avg/max = 54.369/62.456/70.653 ms. `gn_draw_stars` averages 9.059 ms/f; total render averages 9.381 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 2002359 us.

```
isr_wake       1153.4/f 0.263/1.448/15.000 us 2.66%
  isr_pack       144.2/f 5.636/6.775/9.295 us 1.56%
  isr_dma_submit 144.2/f 0.708/0.915/0.980 us 0.21%
```

- Mean pack 6.775 us/call versus submit 0.915 us/call (7.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.66%, corresponding to 60.838 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.167×; Peak live show render: 0.295×; Show end: 0.150×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `gn_draw_stars` — 29.0% of the frame, 18.116 ms/f.
2. `gn_star_scan` — 26.6% of the frame, 16.602 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.140 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared scan-shapes, spherical-math and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
