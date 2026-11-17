# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\mobiusgrid_ship.log`. Replaces `profile_mobiusgrid_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 105 s / 1680 nominal frames at 16 fps; window 16. Raw capture 150 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh MobiusGrid profile 108 16`; analyze only the 105 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:99064, data:152728, headers:8304; free:1771520 / RAM1: variables:314112, code:26008, padding:6760; free:177408 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 689–704 root cycles / 600 MHz = 998171.112 us; measured wall sum 998173 us, delta 1.89 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:48, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `23829075660e13957ccc0f5ace4b87b63212fffb38f9b90d162d47131921bd3f`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 105 s (1680 nominal frames at 16 fps). Show-frame peak render is **55.942 ms**, spilled **0/1680**. Setup frame 1 is **55.942 ms** and remains separate. Live show peak is **33.522 ms**, spilled **0/1679**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1680. The raw log retains 2376 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 33.522 (2), spilled 🟢 0/1679 (0.00%).

README setup frame 1: 55.942 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    59.830 ms  35.898 Mcyc 100.0%
  pov_preserve_half      132.01 us   79.21 kcyc   0.2% x0.9375 141 us/call
  fx_shader_draw         28.538 ms  17.123 Mcyc  47.7%
    fx_prepare_column    690.52 us  414.31 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    94.89 us   56.93 kcyc   0.2% x1 95 us/call
  fx_prepare_frame       818.64 us  491.19 kcyc   1.4% x1 819 us/call
  fx_advance              1.150 ms  690.17 kcyc   1.9%
    palette_cycle_step    1.115 ms  668.88 kcyc   1.9% x1 1115 us/call
  fx_timeline_step        92.54 us   55.52 kcyc   0.2% x1 93 us/call
  canvas_clear            98.71 us   59.22 kcyc   0.2% x1 99 us/call
  canvas_buffer_wait     28.682 ms  17.209 Mcyc  47.9% x1 28682 us/call
```

Wall min/avg/max = 27.847/59.829/66.234 ms. `fx_shader_draw` averages 28.538 ms/f; total render averages 31.148 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 689–704)

```
frame                    62.386 ms  37.431 Mcyc 100.0%
  pov_preserve_half      145.93 us   87.56 kcyc   0.2% x1 146 us/call
  fx_shader_draw         28.523 ms  17.114 Mcyc  45.7%
    fx_prepare_column     1.209 ms  725.67 kcyc   1.9% x144 8 us/call
  fx_palette_base_bake   113.78 us   68.27 kcyc   0.2% x1 114 us/call
  fx_prepare_frame       868.52 us  521.11 kcyc   1.4% x1 869 us/call
  fx_advance              1.011 ms  606.75 kcyc   1.6%
    palette_cycle_step   980.12 us  588.07 kcyc   1.6% x1 980 us/call
  fx_timeline_step       179.36 us  107.62 kcyc   0.3% x1 179 us/call
  canvas_clear            90.00 us   54.00 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     31.247 ms  18.748 Mcyc  50.1% x1 31247 us/call
```

Wall min/avg/max = 58.055/62.385/66.746 ms. `fx_shader_draw` averages 28.523 ms/f; total render averages 31.138 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1665–1680)

```
frame                    62.388 ms  37.433 Mcyc 100.0%
  pov_preserve_half      144.62 us   86.77 kcyc   0.2% x1 145 us/call
  fx_shader_draw         27.130 ms  16.278 Mcyc  43.5%
    fx_prepare_column    775.54 us  465.32 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake   109.38 us   65.63 kcyc   0.2% x1 109 us/call
  fx_prepare_frame       844.83 us  506.90 kcyc   1.4% x1 845 us/call
  fx_advance              1.032 ms  619.00 kcyc   1.7%
    palette_cycle_step   996.80 us  598.08 kcyc   1.6% x1 997 us/call
  fx_timeline_step       100.91 us   60.55 kcyc   0.2% x1 101 us/call
  canvas_clear            91.53 us   54.92 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     32.719 ms  19.631 Mcyc  52.4% x1 32719 us/call
```

Wall min/avg/max = 58.673/62.388/66.033 ms. `fx_shader_draw` averages 27.130 ms/f; total render averages 29.670 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 2 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · Mobius Grid | 37 | 1–16 | held | — | — | 28.538 | 31.148 | 16 | 32.511 | 0/600 |
| 2 · Mobius Grid 2 | 36 | 1377–1392 | held | — | — | 27.148 | 29.630 | 16 | 33.522 | 0/1079 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 689–704; CPU shares use its separately snapshotted ISR interval of 1000219 us.

```
isr_wake       1152.3/f 0.325/1.698/27.486 us 3.13%
  isr_pack       144.0/f 7.343/8.300/10.778 us 1.91%
  isr_dma_submit 144.0/f 0.671/0.913/0.971 us 0.21%
```

- Mean pack 8.300 us/call versus submit 0.913 us/call (9.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.13%, corresponding to 60.544 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.498×; Peak live show render: 0.498×; Show end: 0.475×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 45.7% of the frame, 28.523 ms/f.
2. `fx_prepare_column` — 1.9% of the frame, 1.209 ms/f.
3. `fx_advance` — 1.6% of the frame, 1.011 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MobiusGrid`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
