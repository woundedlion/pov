# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\kaleidoscopeflowers_ship.log`. Replaces `profile_kaleidoscopeflowers_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 173 s / 2768 nominal frames at 16 fps; window 16. Raw capture 218 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh KaleidoscopeFlowers profile 176 16`; analyze only the 173 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:107544, data:153352, headers:8416; free:1762304 / RAM1: variables:314144, code:28888, padding:3880; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 1953–1968 root cycles / 600 MHz = 1002177.295 us; measured wall sum 1002177 us, delta 0.29 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:18, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `366309e33af98cd18ae7e1cac559a076f6de7140eb757cccb8ff5cd6f084e628`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 173 s (2768 nominal frames at 16 fps). Show-frame peak render is **61.934 ms**, spilled **0/2768**. Setup frame 1 is **61.934 ms** and remains separate. Live show peak is **41.644 ms**, spilled **0/2767**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–2768. The raw log retains 3464 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 41.644 (3), spilled 🟢 0/2767 (0.00%).

README setup frame 1: 61.934 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    60.317 ms  36.190 Mcyc 100.0%
  pov_preserve_half      135.67 us   81.40 kcyc   0.2% x0.9375 145 us/call
  fx_shader_draw         29.824 ms  17.895 Mcyc  49.4%
    fx_prepare_column    800.99 us  480.59 kcyc   1.3% x144 6 us/call
  fx_palette_base_bake   772.03 us  463.22 kcyc   1.3%
    fx_palette_rgb_bake  653.08 us  391.85 kcyc   1.1% x1 653 us/call
  fx_prepare_frame       938.79 us  563.27 kcyc   1.6% x1 939 us/call
  fx_advance              1.152 ms  691.34 kcyc   1.9%
    palette_cycle_step    1.102 ms  661.26 kcyc   1.8% x1 1102 us/call
  fx_timeline_step        82.01 us   49.20 kcyc   0.1% x1 82 us/call
  canvas_clear            97.43 us   58.46 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     27.174 ms  16.304 Mcyc  45.1% x1 27174 us/call
```

Wall min/avg/max = 29.809/60.317/66.518 ms. `fx_shader_draw` averages 29.824 ms/f; total render averages 33.143 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 1953–1968)

```
frame                    62.636 ms  37.582 Mcyc 100.0%
  pov_preserve_half      138.48 us   83.09 kcyc   0.2% x1 138 us/call
  fx_shader_draw         33.412 ms  20.047 Mcyc  53.3%
    fx_prepare_column     1.709 ms   1.026 Mcyc   2.7% x144 12 us/call
  fx_palette_base_bake    2.312 ms   1.387 Mcyc   3.7%
    fx_palette_rgb_bake   2.181 ms   1.309 Mcyc   3.5% x1 2181 us/call
  fx_prepare_frame       824.55 us  494.73 kcyc   1.3% x1 825 us/call
  fx_advance             918.13 us  550.88 kcyc   1.5%
    palette_cycle_step   866.48 us  519.89 kcyc   1.4% x1 866 us/call
  fx_timeline_step       181.50 us  108.90 kcyc   0.3% x1 181 us/call
  canvas_clear            89.94 us   53.96 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     24.625 ms  14.775 Mcyc  39.3% x1 24625 us/call
```

Wall min/avg/max = 56.557/62.636/68.060 ms. `fx_shader_draw` averages 33.412 ms/f; total render averages 38.011 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 1633–1648)

```
frame                    62.460 ms  37.476 Mcyc 100.0%
  pov_preserve_half      141.80 us   85.08 kcyc   0.2% x1 142 us/call
  fx_shader_draw         31.951 ms  19.171 Mcyc  51.2%
    fx_prepare_column    803.01 us  481.81 kcyc   1.3% x144 6 us/call
  fx_palette_base_bake    1.452 ms  870.99 kcyc   2.3%
    fx_palette_rgb_bake   1.327 ms  796.08 kcyc   2.1% x1 1327 us/call
  fx_prepare_frame       813.39 us  488.03 kcyc   1.3% x1 813 us/call
  fx_advance             952.20 us  571.32 kcyc   1.5%
    palette_cycle_step   901.05 us  540.63 kcyc   1.4% x1 901 us/call
  fx_timeline_step        81.88 us   49.13 kcyc   0.1% x1 82 us/call
  canvas_clear            91.48 us   54.89 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     26.846 ms  16.107 Mcyc  43.0% x1 26846 us/call
```

Wall min/avg/max = 58.243/62.460/66.583 ms. `fx_shader_draw` averages 31.951 ms/f; total render averages 35.615 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 2753–2768)

```
frame                    62.442 ms  37.465 Mcyc 100.0%
  pov_preserve_half      142.68 us   85.61 kcyc   0.2% x1 143 us/call
  fx_shader_draw         30.068 ms  18.041 Mcyc  48.2%
    fx_prepare_column     1.339 ms  803.24 kcyc   2.1% x144 9 us/call
  fx_palette_base_bake    1.114 ms  668.51 kcyc   1.8%
    fx_palette_rgb_bake  982.68 us  589.61 kcyc   1.6% x1 983 us/call
  fx_prepare_frame       645.98 us  387.59 kcyc   1.0% x1 646 us/call
  fx_advance             949.86 us  569.92 kcyc   1.5%
    palette_cycle_step   890.97 us  534.58 kcyc   1.4% x1 891 us/call
  fx_timeline_step       139.97 us   83.98 kcyc   0.2% x1 140 us/call
  canvas_clear            90.78 us   54.47 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     29.151 ms  17.491 Mcyc  46.7% x1 29151 us/call
```

Wall min/avg/max = 60.051/62.441/65.331 ms. `fx_shader_draw` averages 30.068 ms/f; total render averages 33.291 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 3 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 2 · Open Grid | 36 | 1633–1648 | held | — | — | 31.951 | 35.615 | 16 | 38.595 | 0/1079 |
| 3 · Fine Grid | 37 | 2641–2656 | held | — | — | 31.915 | 35.257 | 16 | 41.644 | 0/1079 |
| 1 · Double Map | 36 | 193–208 | held | — | — | 31.267 | 34.943 | 16 | 38.637 | 0/609 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 1953–1968; CPU shares use its separately snapshotted ISR interval of 1003849 us.

```
isr_wake       1156.5/f 0.288/1.726/23.348 us 3.18%
  isr_pack       144.6/f 7.896/8.593/11.101 us 1.98%
  isr_dma_submit 144.6/f 0.633/0.915/0.985 us 0.21%
```

- Mean pack 8.593 us/call versus submit 0.915 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.18%, corresponding to 60.512 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.530×; Peak live show render: 0.608×; Costliest qualified preset: 0.570×; Show end: 0.533×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 53.3% of the frame, 33.412 ms/f.
2. `fx_palette_base_bake` — 3.7% of the frame, 2.312 ms/f.
3. `fx_palette_rgb_bake` — 3.5% of the frame, 2.181 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeFlowers`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
