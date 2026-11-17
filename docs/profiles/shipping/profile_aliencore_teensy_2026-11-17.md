# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\aliencore_ship.log`. Replaces `profile_aliencore_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | AlienCore 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh AlienCore profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:105336, data:152784, headers:9144; free:1764352 / RAM1: variables:314144, code:28664, padding:4104; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 385–416 root cycles / 600 MHz = 1996691.882 us; measured wall sum 1996695 us, delta 1.56 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:33, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `56998ce8138a94b62b551fd1dee6ffae73c2435bd6188d7c268ca1c785b8d698`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **42.245 ms**, spilled **0/480**. Setup frame 1 is **42.245 ms** and remains separate. Live show peak is **27.530 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 27.530, spilled 🟢 0/479 (0.00%).

README setup frame 1: 42.245 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    60.441 ms  36.265 Mcyc 100.0%
  pov_preserve_half      141.43 us   84.86 kcyc   0.2% x0.96875 146 us/call
  fx_shader_draw         19.257 ms  11.554 Mcyc  31.9%
    fx_prepare_column    921.51 us  552.90 kcyc   1.5% x144 6 us/call
  fx_palette_base_bake   594.17 us  356.50 kcyc   1.0%
    fx_palette_rgb_bake  486.89 us  292.14 kcyc   0.8% x1 487 us/call
  fx_prepare_frame       663.26 us  397.95 kcyc   1.1% x1 663 us/call
  fx_advance             927.86 us  556.72 kcyc   1.5%
    palette_cycle_step   891.56 us  534.93 kcyc   1.5% x1 892 us/call
  fx_timeline_step        79.15 us   47.49 kcyc   0.1% x1 79 us/call
  canvas_clear            89.51 us   53.71 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     38.557 ms  23.134 Mcyc  63.8% x1 38557 us/call
```

Wall min/avg/max = 19.860/60.441/65.593 ms. `fx_shader_draw` averages 19.257 ms/f; total render averages 21.884 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 385–416)

```
frame                    62.397 ms  37.438 Mcyc 100.0%
  pov_preserve_half      150.60 us   90.36 kcyc   0.2% x1 151 us/call
  fx_shader_draw         18.957 ms  11.374 Mcyc  30.4%
    fx_prepare_column    930.77 us  558.46 kcyc   1.5% x144 6 us/call
  fx_palette_base_bake    1.474 ms  884.42 kcyc   2.4%
    fx_palette_rgb_bake   1.357 ms  813.97 kcyc   2.2% x1 1357 us/call
  fx_prepare_frame       723.24 us  433.95 kcyc   1.2% x1 723 us/call
  fx_advance             900.48 us  540.29 kcyc   1.4%
    palette_cycle_step   867.48 us  520.49 kcyc   1.4% x1 867 us/call
  fx_timeline_step        82.87 us   49.72 kcyc   0.1% x1 83 us/call
  canvas_clear            85.32 us   51.19 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     39.901 ms  23.941 Mcyc  63.9% x1 39901 us/call
```

Wall min/avg/max = 55.365/62.396/69.513 ms. `fx_shader_draw` averages 18.957 ms/f; total render averages 22.496 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.397 ms  37.438 Mcyc 100.0%
  pov_preserve_half      149.04 us   89.42 kcyc   0.2% x1 149 us/call
  fx_shader_draw         18.781 ms  11.268 Mcyc  30.1%
    fx_prepare_column    920.42 us  552.25 kcyc   1.5% x144 6 us/call
  fx_palette_base_bake    1.156 ms  693.74 kcyc   1.9%
    fx_palette_rgb_bake   1.040 ms  624.28 kcyc   1.7% x1 1040 us/call
  fx_prepare_frame       702.35 us  421.41 kcyc   1.1% x1 702 us/call
  fx_advance             951.95 us  571.17 kcyc   1.5%
    palette_cycle_step   919.95 us  551.97 kcyc   1.5% x1 920 us/call
  fx_timeline_step        82.90 us   49.74 kcyc   0.1% x1 83 us/call
  canvas_clear            85.42 us   51.25 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     40.364 ms  24.219 Mcyc  64.7% x1 40364 us/call
```

Wall min/avg/max = 55.451/62.396/69.310 ms. `fx_shader_draw` averages 18.781 ms/f; total render averages 22.032 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 385–416; CPU shares use its separately snapshotted ISR interval of 1999897 us.

```
isr_wake       1152.0/f 0.393/1.713/13.941 us 3.15%
  isr_pack       144.0/f 7.771/8.558/11.100 us 1.97%
  isr_dma_submit 144.0/f 0.571/0.915/1.076 us 0.21%
```

- Mean pack 8.558 us/call versus submit 0.915 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.15%, corresponding to 60.531 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.350×; Peak live show render: 0.360×; Show end: 0.353×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 30.4% of the frame, 18.957 ms/f.
2. `fx_palette_base_bake` — 2.4% of the frame, 1.474 ms/f.
3. `fx_palette_rgb_bake` — 2.2% of the frame, 1.357 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienCore`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
