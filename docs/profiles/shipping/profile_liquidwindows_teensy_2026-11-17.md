# LiquidWindows on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\liquidwindows_ship.log`. Replaces `profile_liquidwindows_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | LiquidWindows 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh LiquidWindows profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:98000, data:153152, headers:8944; free:1771520 / RAM1: variables:314144, code:25336, padding:7432; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 257–288 root cycles / 600 MHz = 1996939.047 us; measured wall sum 1996931 us, delta 4.03 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:41, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `4245d1748263aa1fe638772ab5ba8671f4825fca49cdb04e45fd8758087ad28c`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **63.633 ms**, spilled **1/480**. Setup frame 1 is **63.633 ms** and remains separate. Live show peak is **42.732 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 42.732, spilled 🟢 0/479 (0.00%).

README setup frame 1: 63.633 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    61.521 ms  36.912 Mcyc 100.0%
  pov_preserve_half      136.72 us   82.03 kcyc   0.2% x0.96875 141 us/call
  fx_shader_draw         31.895 ms  19.137 Mcyc  51.8%
    fx_prepare_column     1.175 ms  705.00 kcyc   1.9% x144 8 us/call
  fx_palette_base_bake   866.02 us  519.61 kcyc   1.4%
    fx_palette_rgb_bake  721.22 us  432.73 kcyc   1.2% x1 721 us/call
  fx_prepare_frame       819.77 us  491.86 kcyc   1.3% x1 820 us/call
  fx_advance             968.75 us  581.25 kcyc   1.6%
    palette_cycle_step   950.88 us  570.53 kcyc   1.5% x1 951 us/call
  fx_timeline_step        15.73 us    9.44 kcyc   0.0% x1 16 us/call
  canvas_clear            97.13 us   58.28 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     26.551 ms  15.930 Mcyc  43.2% x1 26551 us/call
```

Wall min/avg/max = 32.281/61.520/66.500 ms. `fx_shader_draw` averages 31.895 ms/f; total render averages 34.970 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 257–288)

```
frame                    62.404 ms  37.443 Mcyc 100.0%
  pov_preserve_half      142.33 us   85.40 kcyc   0.2% x1 142 us/call
  fx_shader_draw         32.135 ms  19.281 Mcyc  51.5%
    fx_prepare_column     1.187 ms  712.10 kcyc   1.9% x144 8 us/call
  fx_palette_base_bake    1.879 ms   1.128 Mcyc   3.0%
    fx_palette_rgb_bake   1.734 ms   1.040 Mcyc   2.8% x1 1734 us/call
  fx_prepare_frame       910.28 us  546.17 kcyc   1.5% x1 910 us/call
  fx_advance              1.021 ms  612.80 kcyc   1.6%
    palette_cycle_step    1.003 ms  601.99 kcyc   1.6% x1 1003 us/call
  fx_timeline_step        15.55 us    9.33 kcyc   0.0% x1 16 us/call
  canvas_clear            93.31 us   55.98 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     26.038 ms  15.623 Mcyc  41.7% x1 26038 us/call
```

Wall min/avg/max = 53.765/62.404/71.045 ms. `fx_shader_draw` averages 32.135 ms/f; total render averages 36.367 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.406 ms  37.444 Mcyc 100.0%
  pov_preserve_half      138.81 us   83.28 kcyc   0.2% x1 139 us/call
  fx_shader_draw         30.142 ms  18.085 Mcyc  48.3%
    fx_prepare_column     1.206 ms  723.52 kcyc   1.9% x144 8 us/call
  fx_palette_base_bake    1.019 ms  611.33 kcyc   1.6%
    fx_palette_rgb_bake  876.37 us  525.82 kcyc   1.4% x1 876 us/call
  fx_prepare_frame       651.76 us  391.06 kcyc   1.0% x1 652 us/call
  fx_advance              1.013 ms  608.03 kcyc   1.6%
    palette_cycle_step   996.26 us  597.76 kcyc   1.6% x1 996 us/call
  fx_timeline_step        14.39 us    8.63 kcyc   0.0% x1 14 us/call
  canvas_clear            92.36 us   55.41 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     29.161 ms  17.497 Mcyc  46.7% x1 29161 us/call
```

Wall min/avg/max = 59.009/62.405/65.883 ms. `fx_shader_draw` averages 30.142 ms/f; total render averages 33.245 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 257–288; CPU shares use its separately snapshotted ISR interval of 1999696 us.

```
isr_wake       1151.9/f 0.285/1.688/28.275 us 3.11%
  isr_pack       144.0/f 7.736/8.461/11.065 us 1.94%
  isr_dma_submit 144.0/f 0.620/0.915/4.918 us 0.21%
```

- Mean pack 8.461 us/call versus submit 0.915 us/call (9.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.11%, corresponding to 60.556 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.560×; Peak live show render: 0.582×; Show end: 0.532×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 51.5% of the frame, 32.135 ms/f.
2. `fx_palette_base_bake` — 3.0% of the frame, 1.879 ms/f.
3. `fx_palette_rgb_bake` — 2.8% of the frame, 1.734 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=LiquidWindows`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
