# GridWorms on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\gridworms_ship.log`. Replaces `profile_gridworms_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | GridWorms 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh GridWorms profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:99424, data:152744, headers:8952; free:1770496 / RAM1: variables:314144, code:28696, padding:4072; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 289–320 root cycles / 600 MHz = 1996269.228 us; measured wall sum 1996270 us, delta 0.39 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:35, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `9f912e6313845f4292477fd6b4ff9656424ebc310187a5523d4e273fcf2d12dc`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **54.282 ms**, spilled **0/480**. Setup frame 1 is **54.282 ms** and remains separate. Live show peak is **39.075 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 39.075, spilled 🟢 0/479 (0.00%).

README setup frame 1: 54.282 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    61.043 ms  36.626 Mcyc 100.0%
  pov_preserve_half      137.34 us   82.40 kcyc   0.2% x0.96875 142 us/call
  fx_shader_draw         26.610 ms  15.966 Mcyc  43.6%
    fx_prepare_column    816.62 us  489.97 kcyc   1.3% x144 6 us/call
  fx_palette_base_bake   838.24 us  502.94 kcyc   1.4%
    fx_palette_rgb_bake  715.18 us  429.11 kcyc   1.2% x1 715 us/call
  fx_prepare_frame        1.208 ms  724.52 kcyc   2.0% x1 1208 us/call
  fx_advance             935.77 us  561.46 kcyc   1.5%
    palette_cycle_step   899.47 us  539.68 kcyc   1.5% x1 899 us/call
  fx_timeline_step       101.52 us   60.91 kcyc   0.2% x1 102 us/call
  canvas_clear            90.16 us   54.10 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     30.980 ms  18.588 Mcyc  50.8% x1 30980 us/call
```

Wall min/avg/max = 26.806/61.043/70.017 ms. `fx_shader_draw` averages 26.610 ms/f; total render averages 30.063 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 289–320)

```
frame                    62.383 ms  37.430 Mcyc 100.0%
  pov_preserve_half      141.50 us   84.90 kcyc   0.2% x1 142 us/call
  fx_shader_draw         28.315 ms  16.989 Mcyc  45.4%
    fx_prepare_column    884.26 us  530.55 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake    2.335 ms   1.401 Mcyc   3.7%
    fx_palette_rgb_bake   2.208 ms   1.325 Mcyc   3.5% x1 2208 us/call
  fx_prepare_frame        1.283 ms  769.81 kcyc   2.1% x1 1283 us/call
  fx_advance             893.85 us  536.31 kcyc   1.4%
    palette_cycle_step   858.72 us  515.23 kcyc   1.4% x1 859 us/call
  fx_timeline_step       100.00 us   60.00 kcyc   0.2% x1 100 us/call
  canvas_clear            88.86 us   53.31 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     29.083 ms  17.450 Mcyc  46.6% x1 29083 us/call
```

Wall min/avg/max = 53.391/62.383/70.729 ms. `fx_shader_draw` averages 28.315 ms/f; total render averages 33.300 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.356 ms  37.413 Mcyc 100.0%
  pov_preserve_half      143.85 us   86.31 kcyc   0.2% x1 144 us/call
  fx_shader_draw         24.596 ms  14.758 Mcyc  39.4%
    fx_prepare_column    842.07 us  505.24 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake   987.17 us  592.30 kcyc   1.6%
    fx_palette_rgb_bake  861.31 us  516.78 kcyc   1.4% x1 861 us/call
  fx_prepare_frame       942.69 us  565.62 kcyc   1.5% x1 943 us/call
  fx_advance             940.68 us  564.41 kcyc   1.5%
    palette_cycle_step   904.39 us  542.64 kcyc   1.5% x1 904 us/call
  fx_timeline_step        99.61 us   59.77 kcyc   0.2% x1 100 us/call
  canvas_clear            85.27 us   51.16 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     34.418 ms  20.651 Mcyc  55.2% x1 34418 us/call
```

Wall min/avg/max = 58.842/62.355/65.862 ms. `fx_shader_draw` averages 24.596 ms/f; total render averages 27.938 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 289–320; CPU shares use its separately snapshotted ISR interval of 1999661 us.

```
isr_wake       1151.9/f 0.395/1.681/25.655 us 3.10%
  isr_pack       144.0/f 6.243/8.161/10.921 us 1.88%
  isr_dma_submit 144.0/f 0.708/0.920/6.451 us 0.21%
```

- Mean pack 8.161 us/call versus submit 0.920 us/call (8.9×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.10%, corresponding to 60.562 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.481×; Peak live show render: 0.533×; Show end: 0.447×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 45.4% of the frame, 28.315 ms/f.
2. `fx_palette_base_bake` — 3.7% of the frame, 2.335 ms/f.
3. `fx_palette_rgb_bake` — 3.5% of the frame, 2.208 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridWorms`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
