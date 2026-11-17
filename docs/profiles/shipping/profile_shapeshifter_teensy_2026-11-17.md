# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\shapeshifter_ship.log`. Replaces `profile_shapeshifter_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared Plot raster, column-walk and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 269 s / 4304 nominal frames at 16 fps; window 16. Raw capture 285 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh ShapeShifter profile 272 16`; analyze only the 269 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:196216, data:152392, headers:8768; free:1674240 / RAM1: variables:314048, code:38696, padding:26840; free:144704 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 3553–3568 root cycles / 600 MHz = 996185.123 us; measured wall sum 996187 us, delta 1.88 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:00, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `f1a6d7800f8abf9d81fdc6bc1ee04a55d1166ecaf60697f42b54b4aab786264b`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 269 s (4304 nominal frames at 16 fps). Show-frame peak render is **76.786 ms**, spilled **1/4304**. Setup frame 1 is **76.786 ms** and remains separate. Live show peak is **56.021 ms**, spilled **0/4303**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–4304. The raw log retains 4536 frame records; the 232 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 56.021 (9), spilled 🟢 0/4303 (0.00%).

README setup frame 1: 76.786 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    61.899 ms  37.140 Mcyc 100.0%
  pov_preserve_half      131.67 us   79.00 kcyc   0.2% x0.9375 140 us/call
  ss_draw_all            42.766 ms  25.660 Mcyc  69.1%
    ss_surface_knots      5.178 ms   3.107 Mcyc   8.4%
      ss_surface_mask     3.343 ms   2.006 Mcyc   5.4% x20.25 165 us/call
  ss_timeline_step       128.87 us   77.32 kcyc   0.2% x1 129 us/call
  ss_buffer_wait         18.839 ms  11.304 Mcyc  30.4%
    canvas_clear          98.56 us   59.14 kcyc   0.2% x1 99 us/call
    canvas_buffer_wait   18.708 ms  11.225 Mcyc  30.2% x1 18708 us/call
```

Wall min/avg/max = 42.967/61.899/76.787 ms. `ss_draw_all` averages 42.766 ms/f; total render averages 43.060 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 3553–3568)

```
frame                    62.262 ms  37.357 Mcyc 100.0%
  pov_preserve_half      137.60 us   82.56 kcyc   0.2% x1 138 us/call
  ss_draw_all            53.059 ms  31.835 Mcyc  85.2% x1 53059 us/call
  ss_timeline_step        82.64 us   49.58 kcyc   0.1% x1 83 us/call
  ss_buffer_wait          8.951 ms   5.371 Mcyc  14.4%
    canvas_clear          89.54 us   53.73 kcyc   0.1% x1 90 us/call
    canvas_buffer_wait    8.848 ms   5.309 Mcyc  14.2% x1 8848 us/call
```

Wall min/avg/max = 56.893/62.261/67.404 ms. `ss_draw_all` averages 53.059 ms/f; total render averages 53.310 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 3329–3344)

```
frame                    62.622 ms  37.573 Mcyc 100.0%
  pov_preserve_half      138.68 us   83.21 kcyc   0.2% x1 139 us/call
  ss_draw_all            45.242 ms  27.145 Mcyc  72.2% x1 45242 us/call
  ss_timeline_step        86.36 us   51.82 kcyc   0.1% x1 86 us/call
  ss_buffer_wait         17.121 ms  10.273 Mcyc  27.3%
    canvas_clear          88.16 us   52.89 kcyc   0.1% x1 88 us/call
    canvas_buffer_wait   17.018 ms  10.211 Mcyc  27.2% x1 17018 us/call
```

Wall min/avg/max = 58.983/62.621/66.249 ms. `ss_draw_all` averages 45.242 ms/f; total render averages 45.501 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 4289–4304)

```
frame                    63.346 ms  38.007 Mcyc 100.0%
  pov_preserve_half      141.11 us   84.67 kcyc   0.2% x1 141 us/call
  ss_draw_all            30.085 ms  18.051 Mcyc  47.5%
    ss_surface_knots     529.22 us  317.53 kcyc   0.8%
      ss_surface_mask    341.92 us  205.15 kcyc   0.5% x2 171 us/call
  ss_timeline_step       137.18 us   82.31 kcyc   0.2% x1 137 us/call
  ss_buffer_wait         32.943 ms  19.766 Mcyc  52.0%
    canvas_clear          92.56 us   55.54 kcyc   0.1% x1 93 us/call
    canvas_buffer_wait   32.809 ms  19.686 Mcyc  51.8% x1 32809 us/call
```

Wall min/avg/max = 59.143/63.345/77.373 ms. `ss_draw_all` averages 30.085 ms/f; total render averages 30.402 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 9 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 8 | 28 | 3553–3568 | held | — | — | 53.059 | 53.310 | 16 | 56.021 | 0/479 |
| 7 | 28 | 3137–3152 | held | — | — | 50.269 | 50.518 | 16 | 52.882 | 0/479 |
| 1 | 28 | 81–96 | held | — | — | 47.161 | 47.465 | 8 | 53.189 | 0/471 |
| 6 | 28 | 2721–2736 | held | — | — | 42.630 | 42.887 | 16 | 44.150 | 0/479 |
| 4 | 28 | 1809–1824 | held | — | — | 35.211 | 35.520 | 16 | 37.006 | 0/479 |
| 9 | 28 | 3841–3856 | held | — | — | 29.391 | 29.703 | 16 | 31.193 | 0/479 |
| 5 | 28 | 2193–2208 | held | — | — | 26.548 | 26.832 | 16 | 33.631 | 0/479 |
| 2 | 28 | 625–640 | held | — | — | 24.946 | 25.199 | 16 | 29.329 | 0/479 |
| 3 | 28 | 1217–1232 | held | — | — | 22.046 | 22.346 | 16 | 27.349 | 0/479 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 3553–3568; CPU shares use its separately snapshotted ISR interval of 997971 us.

```
isr_wake       1149.7/f 0.288/1.703/24.853 us 3.14%
  isr_pack       143.8/f 7.281/8.125/10.741 us 1.87%
  isr_dma_submit 143.8/f 0.570/0.911/1.178 us 0.21%
```

- Mean pack 8.125 us/call versus submit 0.911 us/call (8.9×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.14%, corresponding to 60.538 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.689×; Peak live show render: 0.853×; Transition / non-hold regime: 0.728×; Show end: 0.486×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `ss_draw_all` — 85.2% of the frame, 53.059 ms/f.
2. `pov_preserve_half` — 0.2% of the frame, 0.138 ms/f.
3. `canvas_clear` — 0.1% of the frame, 0.090 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared Plot raster, column-walk and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
