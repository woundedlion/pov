# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\millefiori_ship.log`. Replaces `profile_millefiori_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Millefiori 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Millefiori profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:110192, data:153620, headers:8572; free:1759232 / RAM1: variables:314112, code:33288, padding:32248; free:144640 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 385–416 root cycles / 600 MHz = 1998450.765 us; measured wall sum 1998450 us, delta 0.38 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:35, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `24fbb93397e842cfec13d587702d839487492f472d8164706e5d20ce400a6d66`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **79.856 ms**, spilled **1/960**. Setup frame 1 is **79.856 ms** and remains separate. Live show peak is **52.832 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 52.832, spilled 🟢 0/959 (0.00%).

README setup frame 1: 79.856 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.369 ms  37.421 Mcyc 100.0%
  pov_preserve_half      133.46 us   80.08 kcyc   0.2% x0.96875 138 us/call
  fx_shader_draw         42.999 ms  25.799 Mcyc  68.9%
    fx_prepare_column     1.407 ms  844.35 kcyc   2.3% x144 10 us/call
  fx_palette_base_bake   123.50 us   74.10 kcyc   0.2% x1 123 us/call
  fx_prepare_frame        1.243 ms  745.66 kcyc   2.0% x1 1243 us/call
  fx_advance             930.97 us  558.58 kcyc   1.5%
    palette_cycle_step   883.19 us  529.91 kcyc   1.4% x1 883 us/call
  fx_timeline_step        83.17 us   49.90 kcyc   0.1% x1 83 us/call
  canvas_clear            95.27 us   57.16 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     16.625 ms   9.975 Mcyc  26.7% x1 16625 us/call
```

Wall min/avg/max = 40.088/62.368/79.860 ms. `fx_shader_draw` averages 42.999 ms/f; total render averages 45.744 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 385–416)

```
frame                    62.452 ms  37.471 Mcyc 100.0%
  pov_preserve_half      136.59 us   81.96 kcyc   0.2% x1 137 us/call
  fx_shader_draw         44.365 ms  26.619 Mcyc  71.0%
    fx_prepare_column     1.408 ms  844.80 kcyc   2.3% x144 10 us/call
  fx_palette_base_bake   125.39 us   75.23 kcyc   0.2% x1 125 us/call
  fx_prepare_frame        1.067 ms  640.06 kcyc   1.7% x1 1067 us/call
  fx_advance              1.011 ms  606.71 kcyc   1.6%
    palette_cycle_step   960.69 us  576.42 kcyc   1.5% x1 961 us/call
  fx_timeline_step        83.57 us   50.14 kcyc   0.1% x1 84 us/call
  canvas_clear            91.07 us   54.64 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     15.437 ms   9.262 Mcyc  24.7% x1 15437 us/call
```

Wall min/avg/max = 53.159/62.451/71.225 ms. `fx_shader_draw` averages 44.365 ms/f; total render averages 47.014 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.441 ms  37.464 Mcyc 100.0%
  pov_preserve_half      136.80 us   82.08 kcyc   0.2% x1 137 us/call
  fx_shader_draw         46.356 ms  27.814 Mcyc  74.2%
    fx_prepare_column     1.405 ms  842.95 kcyc   2.2% x144 10 us/call
  fx_palette_base_bake   121.91 us   73.14 kcyc   0.2% x1 122 us/call
  fx_prepare_frame        1.313 ms  788.03 kcyc   2.1% x1 1313 us/call
  fx_advance             926.77 us  556.06 kcyc   1.5%
    palette_cycle_step   875.33 us  525.20 kcyc   1.4% x1 875 us/call
  fx_timeline_step        80.88 us   48.53 kcyc   0.1% x1 81 us/call
  canvas_clear            91.45 us   54.87 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     13.282 ms   7.969 Mcyc  21.3% x1 13282 us/call
```

Wall min/avg/max = 57.117/62.440/67.682 ms. `fx_shader_draw` averages 46.356 ms/f; total render averages 49.159 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 385–416; CPU shares use its separately snapshotted ISR interval of 2001869 us.

```
isr_wake       1153.2/f 0.285/1.731/26.715 us 3.19%
  isr_pack       144.2/f 7.851/8.550/11.295 us 1.97%
  isr_dma_submit 144.2/f 0.588/0.918/8.341 us 0.21%
```

- Mean pack 8.550 us/call versus submit 0.918 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.19%, corresponding to 60.506 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.732×; Peak live show render: 0.752×; Show end: 0.787×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 71.0% of the frame, 44.365 ms/f.
2. `fx_prepare_column` — 2.3% of the frame, 1.408 ms/f.
3. `fx_prepare_frame` — 1.7% of the frame, 1.067 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Millefiori`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
