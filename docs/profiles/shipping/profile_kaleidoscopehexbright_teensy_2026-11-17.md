# KaleidoscopeHexBright on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\kaleidoscopehexbright_ship.log`. Replaces `profile_kaleidoscopehexbright_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | KaleidoscopeHexBright 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 105 s / 1680 nominal frames at 16 fps; window 16. Raw capture 150 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh KaleidoscopeHexBright profile 108 16`; analyze only the 105 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:105728, data:152608, headers:8928; free:1764352 / RAM1: variables:314144, code:27976, padding:4792; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 385–400 root cycles / 600 MHz = 1000263.268 us; measured wall sum 1000263 us, delta 0.27 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:29, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `65b95c48ff60403907728f0c667c22d0a017ef32f90b477cf129138ac6a316ff`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 105 s (1680 nominal frames at 16 fps). Show-frame peak render is **57.993 ms**, spilled **0/1680**. Setup frame 1 is **57.993 ms** and remains separate. Live show peak is **41.633 ms**, spilled **0/1679**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1680. The raw log retains 2376 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 41.633 (2), spilled 🟢 0/1679 (0.00%).

README setup frame 1: 57.993 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    60.104 ms  36.062 Mcyc 100.0%
  pov_preserve_half      132.94 us   79.77 kcyc   0.2% x0.9375 142 us/call
  fx_shader_draw         27.142 ms  16.285 Mcyc  45.2%
    fx_prepare_column    712.41 us  427.45 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake   786.49 us  471.89 kcyc   1.3%
    fx_palette_rgb_bake  663.17 us  397.90 kcyc   1.1% x1 663 us/call
  fx_prepare_frame        3.910 ms   2.346 Mcyc   6.5% x1 3910 us/call
  fx_advance             991.41 us  594.85 kcyc   1.6%
    palette_cycle_step   941.85 us  565.11 kcyc   1.6% x1 942 us/call
  fx_timeline_step       102.94 us   61.76 kcyc   0.2% x1 103 us/call
  canvas_clear            99.17 us   59.50 kcyc   0.2% x1 99 us/call
  canvas_buffer_wait     26.773 ms  16.064 Mcyc  44.5% x1 26773 us/call
```

Wall min/avg/max = 31.208/60.103/64.899 ms. `fx_shader_draw` averages 27.142 ms/f; total render averages 33.331 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 385–400)

```
frame                    62.516 ms  37.510 Mcyc 100.0%
  pov_preserve_half      142.96 us   85.77 kcyc   0.2% x1 143 us/call
  fx_shader_draw         29.623 ms  17.774 Mcyc  47.4%
    fx_prepare_column    738.75 us  443.25 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    2.219 ms   1.331 Mcyc   3.5%
    fx_palette_rgb_bake   2.093 ms   1.256 Mcyc   3.3% x1 2093 us/call
  fx_prepare_frame        3.830 ms   2.298 Mcyc   6.1% x1 3830 us/call
  fx_advance             988.40 us  593.04 kcyc   1.6%
    palette_cycle_step   943.88 us  566.33 kcyc   1.5% x1 944 us/call
  fx_timeline_step       107.49 us   64.50 kcyc   0.2% x1 107 us/call
  canvas_clear            90.03 us   54.02 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     25.361 ms  15.217 Mcyc  40.6% x1 25361 us/call
```

Wall min/avg/max = 55.161/62.516/70.117 ms. `fx_shader_draw` averages 29.623 ms/f; total render averages 37.155 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 401–416)

```
frame                    62.339 ms  37.404 Mcyc 100.0%
  pov_preserve_half      140.61 us   84.37 kcyc   0.2% x1 141 us/call
  fx_shader_draw         29.774 ms  17.864 Mcyc  47.8%
    fx_prepare_column    731.39 us  438.83 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    1.892 ms   1.135 Mcyc   3.0%
    fx_palette_rgb_bake   1.764 ms   1.058 Mcyc   2.8% x1 1764 us/call
  fx_prepare_frame        3.802 ms   2.281 Mcyc   6.1% x1 3802 us/call
  fx_advance             991.50 us  594.90 kcyc   1.6%
    palette_cycle_step   942.28 us  565.37 kcyc   1.5% x1 942 us/call
  fx_timeline_step       110.17 us   66.10 kcyc   0.2% x1 110 us/call
  canvas_clear            89.34 us   53.60 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     25.381 ms  15.228 Mcyc  40.7% x1 25381 us/call
```

Wall min/avg/max = 55.956/62.339/68.583 ms. `fx_shader_draw` averages 29.774 ms/f; total render averages 36.959 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 913–928)

```
frame                    62.410 ms  37.446 Mcyc 100.0%
  pov_preserve_half      142.26 us   85.36 kcyc   0.2% x1 142 us/call
  fx_shader_draw         30.007 ms  18.004 Mcyc  48.1%
    fx_prepare_column     1.048 ms  628.52 kcyc   1.7% x144 7 us/call
  fx_palette_base_bake    1.285 ms  771.25 kcyc   2.1%
    fx_palette_rgb_bake   1.153 ms  691.57 kcyc   1.8% x1 1153 us/call
  fx_prepare_frame        3.882 ms   2.329 Mcyc   6.2% x1 3882 us/call
  fx_advance             964.08 us  578.45 kcyc   1.5%
    palette_cycle_step   915.27 us  549.16 kcyc   1.5% x1 915 us/call
  fx_timeline_step       188.24 us  112.94 kcyc   0.3% x1 188 us/call
  canvas_clear            91.37 us   54.82 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     25.692 ms  15.415 Mcyc  41.2% x1 25692 us/call
```

Wall min/avg/max = 56.840/62.410/68.033 ms. `fx_shader_draw` averages 30.007 ms/f; total render averages 36.718 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1665–1680)

```
frame                    62.425 ms  37.455 Mcyc 100.0%
  pov_preserve_half      140.95 us   84.57 kcyc   0.2% x1 141 us/call
  fx_shader_draw         28.495 ms  17.097 Mcyc  45.6%
    fx_prepare_column    738.46 us  443.07 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    1.516 ms  909.68 kcyc   2.4%
    fx_palette_rgb_bake   1.390 ms  833.79 kcyc   2.2% x1 1390 us/call
  fx_prepare_frame        4.001 ms   2.400 Mcyc   6.4% x1 4001 us/call
  fx_advance             944.99 us  566.99 kcyc   1.5%
    palette_cycle_step   888.45 us  533.07 kcyc   1.4% x1 888 us/call
  fx_timeline_step       114.07 us   68.44 kcyc   0.2% x1 114 us/call
  canvas_clear            88.48 us   53.09 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     26.967 ms  16.180 Mcyc  43.2% x1 26967 us/call
```

Wall min/avg/max = 61.025/62.425/63.893 ms. `fx_shader_draw` averages 28.495 ms/f; total render averages 35.459 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 2 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · Hex Twin Wave | 37 | 401–416 | held | — | — | 29.774 | 36.959 | 16 | 41.633 | 0/600 |
| 2 · Hex Twin Wave Alt | 36 | 1345–1360 | held | — | — | 28.848 | 35.249 | 16 | 39.640 | 0/1079 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 385–400; CPU shares use its separately snapshotted ISR interval of 1001953 us.

```
isr_wake       1154.3/f 0.398/1.735/27.080 us 3.19%
  isr_pack       144.3/f 7.706/8.563/10.963 us 1.97%
  isr_dma_submit 144.3/f 0.613/0.915/3.840 us 0.21%
```

- Mean pack 8.563 us/call versus submit 0.915 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.19%, corresponding to 60.506 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.533×; Peak live show render: 0.594×; Costliest qualified preset: 0.591×; Transition / non-hold regime: 0.587×; Show end: 0.567×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 47.4% of the frame, 29.623 ms/f.
2. `fx_prepare_frame` — 6.1% of the frame, 3.830 ms/f.
3. `fx_palette_base_bake` — 3.5% of the frame, 2.219 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexBright`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
