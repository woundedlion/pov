# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\kaleidoscopesmooth_ship.log`. Replaces `profile_kaleidoscopesmooth_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 240 s / 3840 nominal frames at 16 fps; window 16. Raw capture 285 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeSmooth profile 243 16`; analyze only the 240 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:107256, data:153348, headers:8708; free:1762304 / RAM1: variables:314144, code:28888, padding:3880; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 977–992 root cycles / 600 MHz = 998304.927 us; measured wall sum 998302 us, delta 2.93 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:14, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `6f8a56a9526f85af4b6a73d4d4813f06dbec9998e7dac08746349fb571974ce1`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 240 s (3840 nominal frames at 16 fps). Show-frame peak render is **63.892 ms**, spilled **1/3840**. Setup frame 1 is **63.892 ms** and remains separate. Live show peak is **40.288 ms**, spilled **0/3839**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–3840. The raw log retains 4536 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 40.288 (4), spilled 🟢 0/3839 (0.00%).

README setup frame 1: 63.892 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    60.409 ms  36.245 Mcyc 100.0%
  pov_preserve_half      133.97 us   80.38 kcyc   0.2% x0.9375 143 us/call
  fx_shader_draw         30.231 ms  18.138 Mcyc  50.0%
    fx_prepare_column    829.20 us  497.52 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake   778.55 us  467.13 kcyc   1.3%
    fx_palette_rgb_bake  652.01 us  391.20 kcyc   1.1% x1 652 us/call
  fx_prepare_frame        1.011 ms  606.88 kcyc   1.7% x1 1011 us/call
  fx_advance              1.177 ms  706.18 kcyc   1.9%
    palette_cycle_step    1.108 ms  665.09 kcyc   1.8% x1 1108 us/call
  fx_timeline_step        82.54 us   49.52 kcyc   0.1% x1 83 us/call
  canvas_clear            98.64 us   59.19 kcyc   0.2% x1 99 us/call
  canvas_buffer_wait     26.757 ms  16.054 Mcyc  44.3% x1 26757 us/call
```

Wall min/avg/max = 30.966/60.408/65.770 ms. `fx_shader_draw` averages 30.231 ms/f; total render averages 33.651 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 977–992)

```
frame                    62.394 ms  37.436 Mcyc 100.0%
  pov_preserve_half      141.55 us   84.93 kcyc   0.2% x1 142 us/call
  fx_shader_draw         33.388 ms  20.033 Mcyc  53.5%
    fx_prepare_column    974.45 us  584.67 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake    1.262 ms  757.36 kcyc   2.0%
    fx_palette_rgb_bake   1.131 ms  678.82 kcyc   1.8% x1 1131 us/call
  fx_prepare_frame       803.25 us  481.95 kcyc   1.3% x1 803 us/call
  fx_advance             938.08 us  562.85 kcyc   1.5%
    palette_cycle_step   878.80 us  527.28 kcyc   1.4% x1 879 us/call
  fx_timeline_step       170.52 us  102.31 kcyc   0.3% x1 171 us/call
  canvas_clear            90.54 us   54.33 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     25.467 ms  15.280 Mcyc  40.8% x1 25467 us/call
```

Wall min/avg/max = 57.001/62.393/67.998 ms. `fx_shader_draw` averages 33.388 ms/f; total render averages 36.927 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 81–96)

```
frame                    62.391 ms  37.435 Mcyc 100.0%
  pov_preserve_half      138.47 us   83.08 kcyc   0.2% x1 138 us/call
  fx_shader_draw         31.928 ms  19.157 Mcyc  51.2%
    fx_prepare_column    871.93 us  523.16 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake    1.052 ms  631.03 kcyc   1.7%
    fx_palette_rgb_bake  918.78 us  551.27 kcyc   1.5% x1 919 us/call
  fx_prepare_frame       850.91 us  510.54 kcyc   1.4% x1 851 us/call
  fx_advance              1.072 ms  643.19 kcyc   1.7%
    palette_cycle_step    1.002 ms  601.25 kcyc   1.6% x1 1002 us/call
  fx_timeline_step        80.64 us   48.39 kcyc   0.1% x1 81 us/call
  canvas_clear            92.38 us   55.43 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     27.053 ms  16.232 Mcyc  43.4% x1 27053 us/call
```

Wall min/avg/max = 56.217/62.391/68.494 ms. `fx_shader_draw` averages 31.928 ms/f; total render averages 35.339 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 833–848)

```
frame                    62.429 ms  37.458 Mcyc 100.0%
  pov_preserve_half      138.15 us   82.89 kcyc   0.2% x1 138 us/call
  fx_shader_draw         34.295 ms  20.577 Mcyc  54.9%
    fx_prepare_column    968.48 us  581.09 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake   975.06 us  585.04 kcyc   1.6%
    fx_palette_rgb_bake  836.74 us  502.04 kcyc   1.3% x1 837 us/call
  fx_prepare_frame       627.42 us  376.45 kcyc   1.0% x1 627 us/call
  fx_advance             967.62 us  580.57 kcyc   1.5%
    palette_cycle_step   902.88 us  541.73 kcyc   1.4% x1 903 us/call
  fx_timeline_step       172.17 us  103.30 kcyc   0.3% x1 172 us/call
  canvas_clear            90.68 us   54.41 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     25.035 ms  15.021 Mcyc  40.1% x1 25035 us/call
```

Wall min/avg/max = 59.603/62.429/65.749 ms. `fx_shader_draw` averages 34.295 ms/f; total render averages 37.395 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 3825–3840)

```
frame                    62.741 ms  37.644 Mcyc 100.0%
  pov_preserve_half      142.61 us   85.57 kcyc   0.2% x1 143 us/call
  fx_shader_draw         26.908 ms  16.145 Mcyc  42.9%
    fx_prepare_column     1.140 ms  684.00 kcyc   1.8% x144 8 us/call
  fx_palette_base_bake   793.22 us  475.93 kcyc   1.3%
    fx_palette_rgb_bake  659.41 us  395.64 kcyc   1.1% x1 659 us/call
  fx_prepare_frame       538.71 us  323.23 kcyc   0.9% x1 539 us/call
  fx_advance             932.01 us  559.21 kcyc   1.5%
    palette_cycle_step   859.21 us  515.52 kcyc   1.4% x1 859 us/call
  fx_timeline_step       103.12 us   61.87 kcyc   0.2% x1 103 us/call
  canvas_clear            91.63 us   54.98 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     33.104 ms  19.863 Mcyc  52.8% x1 33104 us/call
```

Wall min/avg/max = 59.794/62.740/65.589 ms. `fx_shader_draw` averages 26.908 ms/f; total render averages 29.636 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held windows use fixed authored parameters while the scene keeps animating. Qualified windows exclude incoming parameter blends, fade-in and fade-out, and mixed-owner windows.

The show-scoped data covers 4 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · Coupled Grid | 36 | 81–96 | held | — | — | 31.928 | 35.339 | 16 | 40.256 | 0/602 |
| 2 · Direct Grid | 36 | 1569–1584 | held | — | — | 30.361 | 34.080 | 16 | 40.288 | 0/1079 |
| 3 · Double Map | 37 | 2625–2640 | held | — | — | 29.978 | 33.462 | 16 | 38.763 | 0/1079 |
| 4 · Stretched Grid | 36 | 3569–3584 | held | — | — | 28.309 | 31.496 | 16 | 37.977 | 0/1079 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 977–992; CPU shares use its separately snapshotted ISR interval of 1000237 us.

```
isr_wake       1152.3/f 0.285/1.703/23.335 us 3.14%
  isr_pack       144.0/f 7.716/8.418/11.196 us 1.93%
  isr_dma_submit 144.0/f 0.668/0.916/7.131 us 0.21%
```

- Mean pack 8.418 us/call versus submit 0.916 us/call (9.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.14%, corresponding to 60.538 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.538×; Peak live show render: 0.591×; Costliest qualified preset: 0.565×; Transition / non-hold regime: 0.598×; Show end: 0.474×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 53.5% of the frame, 33.388 ms/f.
2. `fx_palette_base_bake` — 2.0% of the frame, 1.262 ms/f.
3. `fx_palette_rgb_bake` — 1.8% of the frame, 1.131 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
