# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river59-evidence\recover-iri3.log`. Replaces `profile_iridescence_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Iridescence 288×144, single-entry playlist, source `4bf24181947f354011ca64b75c00d722bc4c40d2` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Iridescence profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:109128, data:157208, headers:9120; free:1756160 / RAM1: variables:314144, code:28200, padding:4568; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1092696, data:429268, headers:8900; free:500752 / RAM1: variables:313760, code:185576, padding:11032; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 257–288 root cycles / 600 MHz = 1996387.013 us; measured wall sum 1996389 us, delta 1.00 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 15:15, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `8c0b35c32ee61895894c0c333182da882703916f7cea67c8f826fdec0d35c3cf`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Shared phase shading

The previous shipping renderer peaked at 47.242 ms on COM3. Sharing more code across the show brought that peak to 48.252 ms. The phase stack now executes its pixel shader and signed-byte cube blending in Teensy instruction RAM. A separate run of this placement peaked at 45.513 ms, and this qualification peaked at 45.515 ms with zero missed live frames. The same 480 native quadrants matched the earlier RGB16 capture byte for byte.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **87.094 ms**, spilled **1/480**. Setup frame 1 is **87.094 ms** and remains separate. Live show peak is **45.515 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 520 frame records; the 40 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 45.515, spilled 🟢 0/479 (0.00%).

README setup frame 1: 87.094 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.626 ms  37.575 Mcyc 100.0%
  pov_preserve_half      133.54 us   80.12 kcyc   0.2% x0.96875 138 us/call
  fx_shader_draw         43.453 ms  26.072 Mcyc  69.4%
    fx_prepare_column     1.092 ms  655.32 kcyc   1.7% x144 8 us/call
  fx_palette_base_bake    82.02 us   49.21 kcyc   0.1% x1 82 us/call
  fx_prepare_frame        2.465 ms   1.479 Mcyc   3.9% x1 2465 us/call
  fx_advance              12.99 us    7.79 kcyc   0.0% x1 13 us/call
  fx_timeline_step        15.21 us    9.12 kcyc   0.0% x1 15 us/call
  canvas_clear            91.68 us   55.01 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     16.172 ms   9.703 Mcyc  25.8% x1 16172 us/call
```

Wall min/avg/max = 45.084/62.625/87.098 ms. `fx_shader_draw` averages 43.453 ms/f; total render averages 46.453 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 257–288)

```
frame                    62.387 ms  37.432 Mcyc 100.0%
  pov_preserve_half      138.86 us   83.32 kcyc   0.2% x1 139 us/call
  fx_shader_draw         42.236 ms  25.342 Mcyc  67.7%
    fx_prepare_column     1.086 ms  651.85 kcyc   1.7% x144 8 us/call
  fx_palette_base_bake    81.33 us   48.80 kcyc   0.1% x1 81 us/call
  fx_prepare_frame        2.317 ms   1.390 Mcyc   3.7% x1 2317 us/call
  fx_advance              11.77 us    7.06 kcyc   0.0% x1 12 us/call
  fx_timeline_step        15.35 us    9.21 kcyc   0.0% x1 15 us/call
  canvas_clear            86.79 us   52.07 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     17.305 ms  10.383 Mcyc  27.7% x1 17305 us/call
```

Wall min/avg/max = 61.587/62.387/63.163 ms. `fx_shader_draw` averages 42.236 ms/f; total render averages 45.082 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.388 ms  37.433 Mcyc 100.0%
  pov_preserve_half      136.61 us   81.97 kcyc   0.2% x1 137 us/call
  fx_shader_draw         42.131 ms  25.278 Mcyc  67.5%
    fx_prepare_column     1.092 ms  655.29 kcyc   1.8% x144 8 us/call
  fx_palette_base_bake    82.85 us   49.71 kcyc   0.1% x1 83 us/call
  fx_prepare_frame        2.313 ms   1.388 Mcyc   3.7% x1 2313 us/call
  fx_advance              11.62 us    6.97 kcyc   0.0% x1 12 us/call
  fx_timeline_step        15.20 us    9.12 kcyc   0.0% x1 15 us/call
  canvas_clear            88.02 us   52.81 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     17.414 ms  10.448 Mcyc  27.9% x1 17414 us/call
```

Wall min/avg/max = 61.979/62.387/62.839 ms. `fx_shader_draw` averages 42.131 ms/f; total render averages 44.974 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 257–288; CPU shares use its separately snapshotted ISR interval of 2000002 us.

```
isr_wake       1152.1/f 0.513/1.765/21.040 us 3.25%
  isr_pack       144.0/f 6.828/7.520/10.050 us 1.73%
  isr_dma_submit 144.0/f 0.631/0.915/1.313 us 0.21%
```

- Mean pack 7.520 us/call versus submit 0.915 us/call (8.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.25%, corresponding to 60.469 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.743×; Peak live show render: 0.721×; Show end: 0.720×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 67.7% of the frame, 42.236 ms/f.
2. `fx_prepare_frame` — 3.7% of the frame, 2.317 ms/f.
3. `fx_prepare_column` — 1.7% of the frame, 1.086 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Iridescence`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
