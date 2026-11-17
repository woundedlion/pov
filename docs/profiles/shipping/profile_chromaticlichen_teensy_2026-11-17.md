# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\chromaticlichen_ship.log`. Replaces `profile_chromaticlichen_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh ChromaticLichen profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:106352, data:151924, headers:8988; free:1764352 / RAM1: variables:314144, code:30040, padding:2728; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 225–256 root cycles / 600 MHz = 1996335.717 us; measured wall sum 1996335 us, delta 0.36 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:37, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `14c2d86e957899b395b2663420db67eaca54abd64ce9ddd944b898b7aaac8fbe`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **67.214 ms**, spilled **1/480**. Setup frame 1 is **67.214 ms** and remains separate. Live show peak is **45.730 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 45.730, spilled 🟢 0/479 (0.00%).

README setup frame 1: 67.214 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    61.705 ms  37.023 Mcyc 100.0%
  pov_preserve_half      136.41 us   81.85 kcyc   0.2% x0.96875 141 us/call
  fx_shader_draw         32.279 ms  19.367 Mcyc  52.3%
    fx_prepare_column    967.88 us  580.73 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake   843.85 us  506.31 kcyc   1.4%
    fx_palette_rgb_bake  710.79 us  426.48 kcyc   1.2% x1 711 us/call
  fx_prepare_frame        3.747 ms   2.248 Mcyc   6.1% x1 3747 us/call
  fx_advance             975.18 us  585.11 kcyc   1.6%
    palette_cycle_step   938.04 us  562.83 kcyc   1.5% x1 938 us/call
  fx_timeline_step        80.71 us   48.43 kcyc   0.1% x1 81 us/call
  canvas_clear            96.95 us   58.17 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     23.407 ms  14.044 Mcyc  37.9% x1 23407 us/call
```

Wall min/avg/max = 35.552/61.704/67.215 ms. `fx_shader_draw` averages 32.279 ms/f; total render averages 38.298 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 225–256)

```
frame                    62.385 ms  37.431 Mcyc 100.0%
  pov_preserve_half      140.56 us   84.34 kcyc   0.2% x1 141 us/call
  fx_shader_draw         33.303 ms  19.982 Mcyc  53.4%
    fx_prepare_column    996.20 us  597.72 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake    2.088 ms   1.253 Mcyc   3.3%
    fx_palette_rgb_bake   1.957 ms   1.174 Mcyc   3.1% x1 1957 us/call
  fx_prepare_frame        3.897 ms   2.338 Mcyc   6.2% x1 3897 us/call
  fx_advance             928.25 us  556.95 kcyc   1.5%
    palette_cycle_step   889.95 us  533.97 kcyc   1.4% x1 890 us/call
  fx_timeline_step        79.10 us   47.46 kcyc   0.1% x1 79 us/call
  canvas_clear            91.92 us   55.15 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     21.717 ms  13.030 Mcyc  34.8% x1 21717 us/call
```

Wall min/avg/max = 54.807/62.385/70.053 ms. `fx_shader_draw` averages 33.303 ms/f; total render averages 40.668 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.348 ms  37.409 Mcyc 100.0%
  pov_preserve_half      141.44 us   84.86 kcyc   0.2% x1 141 us/call
  fx_shader_draw         30.805 ms  18.483 Mcyc  49.4%
    fx_prepare_column    982.40 us  589.44 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake   991.36 us  594.82 kcyc   1.6%
    fx_palette_rgb_bake  862.20 us  517.32 kcyc   1.4% x1 862 us/call
  fx_prepare_frame        3.534 ms   2.120 Mcyc   5.7% x1 3534 us/call
  fx_advance             953.24 us  571.94 kcyc   1.5%
    palette_cycle_step   914.09 us  548.45 kcyc   1.5% x1 914 us/call
  fx_timeline_step        78.64 us   47.18 kcyc   0.1% x1 79 us/call
  canvas_clear            92.62 us   55.57 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     25.611 ms  15.366 Mcyc  41.1% x1 25611 us/call
```

Wall min/avg/max = 59.027/62.348/65.719 ms. `fx_shader_draw` averages 30.805 ms/f; total render averages 36.738 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 225–256; CPU shares use its separately snapshotted ISR interval of 2000009 us.

```
isr_wake       1152.1/f 0.355/1.703/20.053 us 3.14%
  isr_pack       144.0/f 7.440/8.260/10.921 us 1.90%
  isr_dma_submit 144.0/f 0.608/0.916/8.533 us 0.21%
```

- Mean pack 8.260 us/call versus submit 0.916 us/call (9.0×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.14%, corresponding to 60.538 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.613×; Peak live show render: 0.651×; Show end: 0.588×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 53.4% of the frame, 33.303 ms/f.
2. `fx_prepare_frame` — 6.2% of the frame, 3.897 ms/f.
3. `fx_palette_base_bake` — 3.3% of the frame, 2.088 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
