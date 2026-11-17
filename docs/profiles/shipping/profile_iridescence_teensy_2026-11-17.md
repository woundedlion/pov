# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\iridescence_ship.log`. Replaces `profile_iridescence_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Iridescence 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh Iridescence profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:112552, data:157268, headers:8708; free:1753088 / RAM1: variables:314144, code:28312, padding:4456; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 1–32 root cycles / 600 MHz = 2006337.863 us; measured wall sum 2006342 us, delta 2.06 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:43, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `43d4f1d666080160b036ab5257360dec5389b6e148a3b460cf78e04fb6d98841`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **90.003 ms**, spilled **1/480**. Setup frame 1 is **90.003 ms** and remains separate. Live show peak is **47.685 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 47.685, spilled 🟢 0/479 (0.00%).

README setup frame 1: 90.003 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.698 ms  37.619 Mcyc 100.0%
  pov_preserve_half      132.55 us   79.53 kcyc   0.2% x0.96875 137 us/call
  fx_shader_draw         43.066 ms  25.840 Mcyc  68.7%
    fx_prepare_column     1.339 ms  803.68 kcyc   2.1% x144 9 us/call
  fx_palette_base_bake    86.40 us   51.84 kcyc   0.1% x1 86 us/call
  fx_prepare_frame        3.204 ms   1.923 Mcyc   5.1% x1 3204 us/call
  fx_advance             935.13 us  561.08 kcyc   1.5%
    palette_cycle_step   910.55 us  546.33 kcyc   1.5% x1 911 us/call
  fx_timeline_step        14.77 us    8.86 kcyc   0.0% x1 15 us/call
  canvas_clear            91.46 us   54.88 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     14.949 ms   8.970 Mcyc  23.8% x1 14949 us/call
```

Wall min/avg/max = 45.329/62.698/90.007 ms. `fx_shader_draw` averages 43.066 ms/f; total render averages 47.749 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.358 ms  37.415 Mcyc 100.0%
  pov_preserve_half      137.40 us   82.44 kcyc   0.2% x1 137 us/call
  fx_shader_draw         41.746 ms  25.048 Mcyc  66.9%
    fx_prepare_column     1.343 ms  805.57 kcyc   2.2% x144 9 us/call
  fx_palette_base_bake    84.72 us   50.83 kcyc   0.1% x1 85 us/call
  fx_prepare_frame        3.023 ms   1.814 Mcyc   4.8% x1 3023 us/call
  fx_advance             901.95 us  541.17 kcyc   1.4%
    palette_cycle_step   877.91 us  526.75 kcyc   1.4% x1 878 us/call
  fx_timeline_step        15.18 us    9.11 kcyc   0.0% x1 15 us/call
  canvas_clear            87.19 us   52.31 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     16.148 ms   9.689 Mcyc  25.9% x1 16148 us/call
```

Wall min/avg/max = 60.451/62.358/64.323 ms. `fx_shader_draw` averages 41.746 ms/f; total render averages 46.211 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 1–32; CPU shares use its separately snapshotted ISR interval of 2419163 us.

```
isr_wake       1393.5/f 0.268/1.450/25.948 us 2.67%
  isr_pack       138.3/f 7.785/8.496/11.226 us 1.55%
  isr_dma_submit 138.3/f 0.565/0.916/5.048 us 0.16%
```

- Mean pack 8.496 us/call versus submit 0.916 us/call (9.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.67%, corresponding to 60.831 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.764×; Show end: 0.739×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 68.7% of the frame, 43.066 ms/f.
2. `fx_prepare_frame` — 5.1% of the frame, 3.204 ms/f.
3. `fx_prepare_column` — 2.1% of the frame, 1.339 ms/f.

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
