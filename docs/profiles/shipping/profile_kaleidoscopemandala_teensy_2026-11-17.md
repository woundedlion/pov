# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\kaleidoscopemandala_ship.log`. Replaces `profile_kaleidoscopemandala_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | KaleidoscopeMandala 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 105 s / 1680 nominal frames at 16 fps; window 16. Raw capture 150 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh KaleidoscopeMandala profile 108 16`; analyze only the 105 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:107432, data:153576, headers:8304; free:1762304 / RAM1: variables:314144, code:28760, padding:4008; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 849–864 root cycles / 600 MHz = 998794.077 us; measured wall sum 998795 us, delta 0.92 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:25, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `fe85d60b0a02f7fb00fc32d9ec6d185c2afb4d9d246d5373daa297705c77a7de`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 105 s (1680 nominal frames at 16 fps). Show-frame peak render is **65.217 ms**, spilled **1/1680**. Setup frame 1 is **65.217 ms** and remains separate. Live show peak is **45.881 ms**, spilled **0/1679**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1680. The raw log retains 2376 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 45.881 (2), spilled 🟢 0/1679 (0.00%).

README setup frame 1: 65.217 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    60.674 ms  36.405 Mcyc 100.0%
  pov_preserve_half      131.60 us   78.96 kcyc   0.2% x0.9375 140 us/call
  fx_shader_draw         32.420 ms  19.452 Mcyc  53.4%
    fx_prepare_column    863.46 us  518.07 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake   683.57 us  410.14 kcyc   1.1%
    fx_palette_rgb_bake  558.50 us  335.10 kcyc   0.9% x1 558 us/call
  fx_prepare_frame       855.93 us  513.56 kcyc   1.4% x1 856 us/call
  fx_advance             958.28 us  574.97 kcyc   1.6%
    palette_cycle_step   921.09 us  552.65 kcyc   1.5% x1 921 us/call
  fx_timeline_step        77.63 us   46.58 kcyc   0.1% x1 78 us/call
  canvas_clear            96.93 us   58.16 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     25.307 ms  15.184 Mcyc  41.7% x1 25307 us/call
```

Wall min/avg/max = 31.835/60.674/65.839 ms. `fx_shader_draw` averages 32.420 ms/f; total render averages 35.367 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 849–864)

```
frame                    62.425 ms  37.455 Mcyc 100.0%
  pov_preserve_half      137.84 us   82.70 kcyc   0.2% x1 138 us/call
  fx_shader_draw         35.083 ms  21.050 Mcyc  56.2%
    fx_prepare_column     1.145 ms  686.87 kcyc   1.8% x144 8 us/call
  fx_palette_base_bake    2.338 ms   1.403 Mcyc   3.7%
    fx_palette_rgb_bake   2.203 ms   1.322 Mcyc   3.5% x1 2203 us/call
  fx_prepare_frame        3.660 ms   2.196 Mcyc   5.9% x1 3660 us/call
  fx_advance             929.82 us  557.89 kcyc   1.5%
    palette_cycle_step   892.13 us  535.28 kcyc   1.4% x1 892 us/call
  fx_timeline_step       177.39 us  106.44 kcyc   0.3% x1 177 us/call
  canvas_clear            88.83 us   53.30 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     19.867 ms  11.920 Mcyc  31.8% x1 19867 us/call
```

Wall min/avg/max = 56.040/62.424/68.858 ms. `fx_shader_draw` averages 35.083 ms/f; total render averages 42.558 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 1649–1664)

```
frame                    62.429 ms  37.457 Mcyc 100.0%
  pov_preserve_half      137.81 us   82.68 kcyc   0.2% x1 138 us/call
  fx_shader_draw         34.384 ms  20.631 Mcyc  55.1%
    fx_prepare_column    866.69 us  520.01 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake    1.223 ms  733.55 kcyc   2.0%
    fx_palette_rgb_bake   1.093 ms  655.98 kcyc   1.8% x1 1093 us/call
  fx_prepare_frame       704.05 us  422.43 kcyc   1.1% x1 704 us/call
  fx_advance             952.90 us  571.74 kcyc   1.5%
    palette_cycle_step   914.83 us  548.90 kcyc   1.5% x1 915 us/call
  fx_timeline_step        70.47 us   42.28 kcyc   0.1% x1 70 us/call
  canvas_clear            87.43 us   52.46 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     24.726 ms  14.836 Mcyc  39.6% x1 24726 us/call
```

Wall min/avg/max = 59.178/62.428/65.547 ms. `fx_shader_draw` averages 34.384 ms/f; total render averages 37.703 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1665–1680)

```
frame                    62.604 ms  37.563 Mcyc 100.0%
  pov_preserve_half      139.13 us   83.48 kcyc   0.2% x1 139 us/call
  fx_shader_draw         34.284 ms  20.570 Mcyc  54.8%
    fx_prepare_column    885.94 us  531.56 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake    1.203 ms  721.73 kcyc   1.9%
    fx_palette_rgb_bake   1.075 ms  644.96 kcyc   1.7% x1 1075 us/call
  fx_prepare_frame       900.35 us  540.21 kcyc   1.4% x1 900 us/call
  fx_advance             960.32 us  576.19 kcyc   1.5%
    palette_cycle_step   913.61 us  548.16 kcyc   1.5% x1 914 us/call
  fx_timeline_step        81.81 us   49.09 kcyc   0.1% x1 82 us/call
  canvas_clear            89.89 us   53.94 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     24.806 ms  14.884 Mcyc  39.6% x1 24806 us/call
```

Wall min/avg/max = 59.387/62.604/65.257 ms. `fx_shader_draw` averages 34.284 ms/f; total render averages 37.798 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 2 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 2 · Cup Hue | 36 | 1649–1664 | held | — | — | 34.384 | 37.703 | 16 | 45.881 | 0/1079 |
| 1 · Wave Mirror | 37 | 545–560 | held | — | — | 34.372 | 37.754 | 8 | 43.329 | 0/600 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 849–864; CPU shares use its separately snapshotted ISR interval of 1000608 us.

```
isr_wake       1152.8/f 0.395/1.751/22.415 us 3.23%
  isr_pack       144.1/f 8.031/8.600/11.078 us 1.98%
  isr_dma_submit 144.1/f 0.606/0.923/9.273 us 0.21%
```

- Mean pack 8.600 us/call versus submit 0.923 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.23%, corresponding to 60.481 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.566×; Peak live show render: 0.681×; Costliest qualified preset: 0.603×; Show end: 0.605×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 56.2% of the frame, 35.083 ms/f.
2. `fx_prepare_frame` — 5.9% of the frame, 3.660 ms/f.
3. `fx_palette_base_bake` — 3.7% of the frame, 2.338 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeMandala`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
