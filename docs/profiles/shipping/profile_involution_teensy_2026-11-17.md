# Involution on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\involution_ship.log`. Replaces `profile_involution_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Involution 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Involution profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:114376, data:159696, headers:8552; free:1748992 / RAM1: variables:314144, code:29496, padding:3272; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 1–32 root cycles / 600 MHz = 1998194.697 us; measured wall sum 1998195 us, delta 0.15 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:39, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `3c685f843474bea2f7300820436966c7b4c9363719b6e00ddf15e80e22101e12`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **82.497 ms**, spilled **1/480**. Setup frame 1 is **82.497 ms** and remains separate. Live show peak is **44.738 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 44.738, spilled 🟢 0/479 (0.00%).

README setup frame 1: 82.497 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.444 ms  37.466 Mcyc 100.0%
  pov_preserve_half      134.12 us   80.47 kcyc   0.2% x0.96875 138 us/call
  fx_shader_draw         43.134 ms  25.880 Mcyc  69.1%
    fx_prepare_column     2.024 ms   1.214 Mcyc   3.2% x144 14 us/call
  fx_palette_base_bake    89.73 us   53.84 kcyc   0.1% x1 90 us/call
  fx_prepare_frame       131.64 us   78.98 kcyc   0.2% x1 132 us/call
  fx_advance             927.94 us  556.77 kcyc   1.5%
    palette_cycle_step   871.74 us  523.04 kcyc   1.4% x1 872 us/call
  fx_timeline_step        53.56 us   32.14 kcyc   0.1% x1 54 us/call
  canvas_clear            96.88 us   58.13 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     17.689 ms  10.614 Mcyc  28.3% x1 17689 us/call
```

Wall min/avg/max = 41.789/62.443/82.501 ms. `fx_shader_draw` averages 43.134 ms/f; total render averages 44.754 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.412 ms  37.447 Mcyc 100.0%
  pov_preserve_half      136.51 us   81.91 kcyc   0.2% x1 137 us/call
  fx_shader_draw         41.755 ms  25.053 Mcyc  66.9%
    fx_prepare_column     2.084 ms   1.250 Mcyc   3.3% x144 14 us/call
  fx_palette_base_bake    87.15 us   52.29 kcyc   0.1% x1 87 us/call
  fx_prepare_frame        39.47 us   23.68 kcyc   0.1% x1 39 us/call
  fx_advance             875.93 us  525.56 kcyc   1.4%
    palette_cycle_step   825.18 us  495.11 kcyc   1.3% x1 825 us/call
  fx_timeline_step        49.91 us   29.95 kcyc   0.1% x1 50 us/call
  canvas_clear            92.02 us   55.21 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     19.198 ms  11.519 Mcyc  30.8% x1 19198 us/call
```

Wall min/avg/max = 60.535/62.412/64.300 ms. `fx_shader_draw` averages 41.755 ms/f; total render averages 43.214 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 1–32; CPU shares use its separately snapshotted ISR interval of 2416725 us.

```
isr_wake       1392.1/f 0.270/1.446/18.826 us 2.66%
  isr_pack       138.1/f 6.798/8.561/11.225 us 1.56%
  isr_dma_submit 138.1/f 0.595/0.915/9.916 us 0.16%
```

- Mean pack 8.561 us/call versus submit 0.915 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.66%, corresponding to 60.838 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.716×; Show end: 0.691×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 69.1% of the frame, 43.134 ms/f.
2. `fx_prepare_column` — 3.2% of the frame, 2.024 ms/f.
3. `fx_advance` — 1.5% of the frame, 0.928 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Involution`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
