# GridWarp on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\gridwarp_ship.log`. Replaces `profile_gridwarp_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | GridWarp 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh GridWarp profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:111096, data:154948, headers:8388; free:1757184 / RAM1: variables:314144, code:29240, padding:3528; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 257–288 root cycles / 600 MHz = 1995729.408 us; measured wall sum 1995728 us, delta 0.71 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:41, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `da08da8127df27c87c1cf5b1bf934114e61f1817b966435c9715a32a27245cc9`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **80.625 ms**, spilled **1/480**. Setup frame 1 is **80.625 ms** and remains separate. Live show peak is **56.343 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 56.343, spilled 🟢 0/479 (0.00%).

README setup frame 1: 80.625 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.415 ms  37.449 Mcyc 100.0%
  pov_preserve_half      134.79 us   80.88 kcyc   0.2% x0.96875 139 us/call
  fx_shader_draw         42.873 ms  25.724 Mcyc  68.7%
    fx_prepare_column     5.052 ms   3.031 Mcyc   8.1% x144 35 us/call
  fx_palette_base_bake   863.63 us  518.18 kcyc   1.4%
    fx_palette_rgb_bake  713.86 us  428.32 kcyc   1.1% x1 714 us/call
  fx_prepare_frame        3.688 ms   2.213 Mcyc   5.9% x1 3688 us/call
  fx_advance              1.070 ms  642.13 kcyc   1.7%
    palette_cycle_step    1.029 ms  617.51 kcyc   1.6% x1 1029 us/call
  fx_timeline_step       105.23 us   63.14 kcyc   0.2% x1 105 us/call
  canvas_clear            99.32 us   59.59 kcyc   0.2% x1 99 us/call
  canvas_buffer_wait     13.438 ms   8.063 Mcyc  21.5% x1 13438 us/call
```

Wall min/avg/max = 44.755/62.415/80.628 ms. `fx_shader_draw` averages 42.873 ms/f; total render averages 48.977 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 257–288)

```
frame                    62.367 ms  37.420 Mcyc 100.0%
  pov_preserve_half      138.67 us   83.20 kcyc   0.2% x1 139 us/call
  fx_shader_draw         42.682 ms  25.609 Mcyc  68.4%
    fx_prepare_column     4.666 ms   2.799 Mcyc   7.5% x144 32 us/call
  fx_palette_base_bake    1.856 ms   1.114 Mcyc   3.0%
    fx_palette_rgb_bake   1.710 ms   1.026 Mcyc   2.7% x1 1710 us/call
  fx_prepare_frame        3.801 ms   2.281 Mcyc   6.1% x1 3801 us/call
  fx_advance             897.52 us  538.51 kcyc   1.4%
    palette_cycle_step   855.27 us  513.16 kcyc   1.4% x1 855 us/call
  fx_timeline_step       113.13 us   67.88 kcyc   0.2% x1 113 us/call
  canvas_clear            95.61 us   57.37 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait     12.645 ms   7.587 Mcyc  20.3% x1 12645 us/call
```

Wall min/avg/max = 53.489/62.366/71.216 ms. `fx_shader_draw` averages 42.682 ms/f; total render averages 49.722 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.405 ms  37.443 Mcyc 100.0%
  pov_preserve_half      140.50 us   84.30 kcyc   0.2% x1 141 us/call
  fx_shader_draw         38.705 ms  23.223 Mcyc  62.0%
    fx_prepare_column     4.180 ms   2.508 Mcyc   6.7% x144 29 us/call
  fx_palette_base_bake    1.010 ms  605.77 kcyc   1.6%
    fx_palette_rgb_bake  862.86 us  517.71 kcyc   1.4% x1 863 us/call
  fx_prepare_frame        3.530 ms   2.118 Mcyc   5.7% x1 3530 us/call
  fx_advance             981.83 us  589.10 kcyc   1.6%
    palette_cycle_step   939.83 us  563.90 kcyc   1.5% x1 940 us/call
  fx_timeline_step       108.74 us   65.24 kcyc   0.2% x1 109 us/call
  canvas_clear            95.93 us   57.56 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait     17.694 ms  10.617 Mcyc  28.4% x1 17694 us/call
```

Wall min/avg/max = 56.962/62.404/68.014 ms. `fx_shader_draw` averages 38.705 ms/f; total render averages 44.711 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 257–288; CPU shares use its separately snapshotted ISR interval of 1999114 us.

```
isr_wake       1151.6/f 0.285/1.735/19.686 us 3.19%
  isr_pack       143.9/f 7.738/8.460/11.053 us 1.94%
  isr_dma_submit 143.9/f 0.585/0.915/8.130 us 0.21%
```

- Mean pack 8.460 us/call versus submit 0.915 us/call (9.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.19%, corresponding to 60.506 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.784×; Peak live show render: 0.796×; Show end: 0.715×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 68.4% of the frame, 42.682 ms/f.
2. `fx_prepare_column` — 7.5% of the frame, 4.666 ms/f.
3. `fx_prepare_frame` — 6.1% of the frame, 3.801 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridWarp`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
