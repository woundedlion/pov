# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\mermaidskin_ship.log`. Replaces `profile_mermaidskin_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | MermaidSkin 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh MermaidSkin profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:106480, data:151920, headers:8864; free:1764352 / RAM1: variables:314144, code:30040, padding:2728; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 225–256 root cycles / 600 MHz = 1996197.233 us; measured wall sum 1996201 us, delta 1.89 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:37, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `dbd8fbfc2052c0085dacce648d1f96f1938ffd287b79dd94b11c5c9b2e3072f6`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **68.306 ms**, spilled **1/480**. Setup frame 1 is **68.306 ms** and remains separate. Live show peak is **45.898 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 45.898, spilled 🟢 0/479 (0.00%).

README setup frame 1: 68.306 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    61.742 ms  37.045 Mcyc 100.0%
  pov_preserve_half      136.75 us   82.05 kcyc   0.2% x0.96875 141 us/call
  fx_shader_draw         32.743 ms  19.646 Mcyc  53.0%
    fx_prepare_column    792.72 us  475.63 kcyc   1.3% x144 6 us/call
  fx_palette_base_bake   816.13 us  489.68 kcyc   1.3%
    fx_palette_rgb_bake  694.09 us  416.46 kcyc   1.1% x1 694 us/call
  fx_prepare_frame        3.768 ms   2.261 Mcyc   6.1% x1 3768 us/call
  fx_advance             991.73 us  595.04 kcyc   1.6%
    palette_cycle_step   956.47 us  573.88 kcyc   1.5% x1 956 us/call
  fx_timeline_step        81.82 us   49.09 kcyc   0.1% x1 82 us/call
  canvas_clear            95.53 us   57.32 kcyc   0.2% x1 96 us/call
  canvas_buffer_wait     22.974 ms  13.784 Mcyc  37.2% x1 22974 us/call
```

Wall min/avg/max = 35.909/61.742/68.307 ms. `fx_shader_draw` averages 32.743 ms/f; total render averages 38.769 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 225–256)

```
frame                    62.381 ms  37.429 Mcyc 100.0%
  pov_preserve_half      140.61 us   84.37 kcyc   0.2% x1 141 us/call
  fx_shader_draw         33.500 ms  20.100 Mcyc  53.7%
    fx_prepare_column    806.26 us  483.76 kcyc   1.3% x144 6 us/call
  fx_palette_base_bake    2.092 ms   1.255 Mcyc   3.4%
    fx_palette_rgb_bake   1.965 ms   1.179 Mcyc   3.1% x1 1965 us/call
  fx_prepare_frame        3.884 ms   2.331 Mcyc   6.2% x1 3884 us/call
  fx_advance             926.04 us  555.62 kcyc   1.5%
    palette_cycle_step   890.49 us  534.30 kcyc   1.4% x1 890 us/call
  fx_timeline_step        80.81 us   48.49 kcyc   0.1% x1 81 us/call
  canvas_clear            91.18 us   54.71 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     21.536 ms  12.921 Mcyc  34.5% x1 21536 us/call
```

Wall min/avg/max = 54.796/62.381/69.940 ms. `fx_shader_draw` averages 33.500 ms/f; total render averages 40.846 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.364 ms  37.419 Mcyc 100.0%
  pov_preserve_half      142.28 us   85.37 kcyc   0.2% x1 142 us/call
  fx_shader_draw         31.084 ms  18.650 Mcyc  49.8%
    fx_prepare_column    810.12 us  486.07 kcyc   1.3% x144 6 us/call
  fx_palette_base_bake   981.69 us  589.02 kcyc   1.6%
    fx_palette_rgb_bake  858.66 us  515.20 kcyc   1.4% x1 859 us/call
  fx_prepare_frame        3.605 ms   2.163 Mcyc   5.8% x1 3605 us/call
  fx_advance             962.31 us  577.39 kcyc   1.5%
    palette_cycle_step   927.93 us  556.76 kcyc   1.5% x1 928 us/call
  fx_timeline_step        79.39 us   47.63 kcyc   0.1% x1 79 us/call
  canvas_clear            90.59 us   54.35 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     25.293 ms  15.176 Mcyc  40.6% x1 25293 us/call
```

Wall min/avg/max = 59.080/62.364/65.814 ms. `fx_shader_draw` averages 31.084 ms/f; total render averages 37.071 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 225–256; CPU shares use its separately snapshotted ISR interval of 1999981 us.

```
isr_wake       1152.1/f 0.386/1.760/21.706 us 3.24%
  isr_pack       144.0/f 7.905/8.686/11.303 us 2.00%
  isr_dma_submit 144.0/f 0.601/0.920/8.788 us 0.21%
```

- Mean pack 8.686 us/call versus submit 0.920 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.24%, corresponding to 60.475 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.620×; Peak live show render: 0.654×; Show end: 0.593×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 53.7% of the frame, 33.500 ms/f.
2. `fx_prepare_frame` — 6.2% of the frame, 3.884 ms/f.
3. `fx_palette_base_bake` — 3.4% of the frame, 2.092 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MermaidSkin`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
