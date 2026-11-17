# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\cosmiceyeball_ship.log`. Replaces `profile_cosmiceyeball_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | CosmicEyeball 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh CosmicEyeball profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:96976, data:152268, headers:8804; free:1773568 / RAM1: variables:314144, code:27304, padding:5464; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 705–736 root cycles / 600 MHz = 1996788.925 us; measured wall sum 1996790 us, delta 0.54 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:43, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `956af1a718ff7e806beae859e0302e940803faf5d7d6c6f3ceb5d028c398f1b8`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **48.727 ms**, spilled **0/960**. Setup frame 1 is **48.727 ms** and remains separate. Live show peak is **30.734 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1097 frame records; the 137 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 30.734, spilled 🟢 0/959 (0.00%).

README setup frame 1: 48.727 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    60.797 ms  36.478 Mcyc 100.0%
  pov_preserve_half      139.69 us   83.81 kcyc   0.2% x0.96875 144 us/call
  fx_shader_draw         23.884 ms  14.330 Mcyc  39.3%
    fx_prepare_column     1.519 ms  911.51 kcyc   2.5% x144 11 us/call
  fx_palette_base_bake   532.16 us  319.30 kcyc   0.9%
    fx_palette_rgb_bake  437.46 us  262.48 kcyc   0.7% x1 437 us/call
  fx_prepare_frame       515.24 us  309.15 kcyc   0.8% x1 515 us/call
  fx_advance              1.895 ms   1.137 Mcyc   3.1%
    palette_cycle_step    1.878 ms   1.127 Mcyc   3.1% x1 1878 us/call
  fx_timeline_step       108.96 us   65.38 kcyc   0.2% x1 109 us/call
  canvas_clear            99.97 us   59.98 kcyc   0.2% x1 100 us/call
  canvas_buffer_wait     33.448 ms  20.069 Mcyc  55.0% x1 33448 us/call
```

Wall min/avg/max = 24.780/60.796/66.909 ms. `fx_shader_draw` averages 23.884 ms/f; total render averages 27.349 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 705–736)

```
frame                    62.400 ms  37.440 Mcyc 100.0%
  pov_preserve_half      143.45 us   86.07 kcyc   0.2% x1 143 us/call
  fx_shader_draw         23.423 ms  14.054 Mcyc  37.5%
    fx_prepare_column     1.518 ms  910.68 kcyc   2.4% x144 11 us/call
  fx_palette_base_bake   707.53 us  424.52 kcyc   1.1%
    fx_palette_rgb_bake  580.28 us  348.17 kcyc   0.9% x1 580 us/call
  fx_prepare_frame       513.17 us  307.90 kcyc   0.8% x1 513 us/call
  fx_advance              2.154 ms   1.292 Mcyc   3.5%
    palette_cycle_step    2.138 ms   1.283 Mcyc   3.4% x1 2138 us/call
  fx_timeline_step       113.03 us   67.82 kcyc   0.2% x1 113 us/call
  canvas_clear            96.14 us   57.69 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait     35.078 ms  21.047 Mcyc  56.2% x1 35078 us/call
```

Wall min/avg/max = 55.919/62.399/68.965 ms. `fx_shader_draw` averages 23.423 ms/f; total render averages 27.321 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.383 ms  37.430 Mcyc 100.0%
  pov_preserve_half      145.33 us   87.20 kcyc   0.2% x1 145 us/call
  fx_shader_draw         23.246 ms  13.948 Mcyc  37.3%
    fx_prepare_column     1.513 ms  907.64 kcyc   2.4% x144 11 us/call
  fx_palette_base_bake   825.12 us  495.07 kcyc   1.3%
    fx_palette_rgb_bake  726.98 us  436.19 kcyc   1.2% x1 727 us/call
  fx_prepare_frame       521.52 us  312.91 kcyc   0.8% x1 522 us/call
  fx_advance              2.000 ms   1.200 Mcyc   3.2%
    palette_cycle_step    1.984 ms   1.190 Mcyc   3.2% x1 1984 us/call
  fx_timeline_step       112.76 us   67.66 kcyc   0.2% x1 113 us/call
  canvas_clear            96.97 us   58.18 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     35.263 ms  21.158 Mcyc  56.5% x1 35263 us/call
```

Wall min/avg/max = 56.750/62.382/68.012 ms. `fx_shader_draw` averages 23.246 ms/f; total render averages 27.120 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 705–736; CPU shares use its separately snapshotted ISR interval of 2000304 us.

```
isr_wake       1152.2/f 0.285/1.705/24.033 us 3.14%
  isr_pack       144.0/f 7.945/8.703/11.798 us 2.00%
  isr_dma_submit 144.0/f 0.676/0.913/5.035 us 0.21%
```

- Mean pack 8.703 us/call versus submit 0.913 us/call (9.5×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.14%, corresponding to 60.538 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.438×; Peak live show render: 0.437×; Show end: 0.434×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 37.5% of the frame, 23.423 ms/f.
2. `fx_advance` — 3.5% of the frame, 2.154 ms/f.
3. `palette_cycle_step` — 3.4% of the frame, 2.138 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=CosmicEyeball`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
