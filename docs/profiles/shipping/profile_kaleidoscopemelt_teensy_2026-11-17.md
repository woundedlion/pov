# KaleidoscopeMelt on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\kaleidoscopemelt_ship.log`. Replaces `profile_kaleidoscopemelt_teensy_2026-11-13.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | KaleidoscopeMelt 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeMelt profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:103176, data:152244, headers:8772; free:1767424 / RAM1: variables:314112, code:30440, padding:2328; free:177408 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 161–192 root cycles / 600 MHz = 1999399.085 us; measured wall sum 1999400 us, delta 0.46 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:45, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `cc931dbb9d2ce53f20a8d4d86c3a3a94d70b84f18eb45bf48498b706c03768ca`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **78.699 ms**, spilled **1/960**. Setup frame 1 is **78.699 ms** and remains separate. Live show peak is **50.074 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 50.074, spilled 🟢 0/959 (0.00%).

README setup frame 1: 78.699 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.315 ms  37.389 Mcyc 100.0%
  pov_preserve_half      135.97 us   81.58 kcyc   0.2% x0.96875 140 us/call
  fx_shader_draw         40.926 ms  24.556 Mcyc  65.7%
    fx_prepare_column    681.81 us  409.09 kcyc   1.1% x144 5 us/call
  fx_palette_base_bake   325.08 us  195.05 kcyc   0.5% x1 325 us/call
  fx_prepare_frame        2.943 ms   1.766 Mcyc   4.7% x1 2943 us/call
  fx_advance              1.161 ms  696.90 kcyc   1.9%
    palette_cycle_step    1.133 ms  679.98 kcyc   1.8% x1 1133 us/call
  fx_timeline_step        72.39 us   43.43 kcyc   0.1% x1 72 us/call
  canvas_clear            99.82 us   59.89 kcyc   0.2% x1 100 us/call
  canvas_buffer_wait     16.530 ms   9.918 Mcyc  26.5% x1 16530 us/call
```

Wall min/avg/max = 43.454/62.314/78.700 ms. `fx_shader_draw` averages 40.926 ms/f; total render averages 45.785 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 161–192)

```
frame                    62.481 ms  37.489 Mcyc 100.0%
  pov_preserve_half      137.91 us   82.74 kcyc   0.2% x1 138 us/call
  fx_shader_draw         43.049 ms  25.829 Mcyc  68.9%
    fx_prepare_column    681.44 us  408.86 kcyc   1.1% x144 5 us/call
  fx_palette_base_bake   337.22 us  202.33 kcyc   0.5% x1 337 us/call
  fx_prepare_frame        2.955 ms   1.773 Mcyc   4.7% x1 2955 us/call
  fx_advance             989.28 us  593.57 kcyc   1.6%
    palette_cycle_step   962.47 us  577.48 kcyc   1.5% x1 962 us/call
  fx_timeline_step        74.92 us   44.95 kcyc   0.1% x1 75 us/call
  canvas_clear            95.42 us   57.25 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     14.718 ms   8.831 Mcyc  23.6% x1 14718 us/call
```

Wall min/avg/max = 60.026/62.481/64.922 ms. `fx_shader_draw` averages 43.049 ms/f; total render averages 47.763 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.383 ms  37.430 Mcyc 100.0%
  pov_preserve_half      138.12 us   82.87 kcyc   0.2% x1 138 us/call
  fx_shader_draw         42.125 ms  25.275 Mcyc  67.5%
    fx_prepare_column    678.35 us  407.01 kcyc   1.1% x144 5 us/call
  fx_palette_base_bake   339.82 us  203.89 kcyc   0.5% x1 340 us/call
  fx_prepare_frame        3.006 ms   1.804 Mcyc   4.8% x1 3006 us/call
  fx_advance              1.043 ms  625.74 kcyc   1.7%
    palette_cycle_step    1.016 ms  609.72 kcyc   1.6% x1 1016 us/call
  fx_timeline_step        73.48 us   44.09 kcyc   0.1% x1 73 us/call
  canvas_clear            97.75 us   58.65 kcyc   0.2% x1 98 us/call
  canvas_buffer_wait     15.439 ms   9.264 Mcyc  24.7% x1 15439 us/call
```

Wall min/avg/max = 58.130/62.383/66.786 ms. `fx_shader_draw` averages 42.125 ms/f; total render averages 46.944 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 161–192; CPU shares use its separately snapshotted ISR interval of 2002246 us.

```
isr_wake       1153.4/f 0.285/1.751/21.681 us 3.23%
  isr_pack       144.2/f 7.911/8.671/11.735 us 1.99%
  isr_dma_submit 144.2/f 0.651/0.916/1.173 us 0.21%
```

- Mean pack 8.671 us/call versus submit 0.916 us/call (9.5×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.23%, corresponding to 60.481 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.733×; Peak live show render: 0.764×; Show end: 0.751×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 68.9% of the frame, 43.049 ms/f.
2. `fx_prepare_frame` — 4.7% of the frame, 2.955 ms/f.
3. `fx_advance` — 1.6% of the frame, 0.989 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeMelt`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
