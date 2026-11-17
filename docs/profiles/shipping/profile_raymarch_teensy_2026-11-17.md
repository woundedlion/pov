# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\raymarch_ship.log`. Replaces `profile_raymarch_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Raymarch profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:108552, data:186576, headers:9000; free:1727488 / RAM1: variables:313920, code:25736, padding:7032; free:177600 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 1–32 root cycles / 600 MHz = 2010778.698 us; measured wall sum 2010776 us, delta 1.34 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:52, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `7576f18dc8230c1575155827c746f53f91620a84fe5ea0c719270208733720a2`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **94.040 ms**, spilled **1/480**. Setup frame 1 is **94.040 ms** and remains separate. Live show peak is **52.759 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 52.759, spilled 🟢 0/479 (0.00%).

README setup frame 1: 94.040 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.837 ms  37.702 Mcyc 100.0%
  pov_preserve_half      130.38 us   78.23 kcyc   0.2% x0.96875 135 us/call
  rm_shader_draw         45.952 ms  27.571 Mcyc  73.1% x1 45952 us/call
  rm_timeline_step       330.53 us  198.32 kcyc   0.5% x1 331 us/call
  canvas_clear            89.95 us   53.97 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     13.419 ms   8.051 Mcyc  21.4% x1 13419 us/call
```

Wall min/avg/max = 51.456/62.836/94.041 ms. `rm_shader_draw` averages 45.952 ms/f; total render averages 49.418 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.277 ms  37.366 Mcyc 100.0%
  pov_preserve_half      133.80 us   80.28 kcyc   0.2% x1 134 us/call
  rm_shader_draw         43.302 ms  25.981 Mcyc  69.5% x1 43302 us/call
  rm_timeline_step       331.15 us  198.69 kcyc   0.5% x1 331 us/call
  canvas_clear            85.40 us   51.24 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     15.499 ms   9.299 Mcyc  24.9% x1 15499 us/call
```

Wall min/avg/max = 56.684/62.277/68.985 ms. `rm_shader_draw` averages 43.302 ms/f; total render averages 46.778 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 1–32; CPU shares use its separately snapshotted ISR interval of 2417721 us.

```
isr_wake       1392.7/f 0.253/1.321/29.828 us 2.43%
  isr_pack       138.2/f 5.511/7.236/10.153 us 1.32%
  isr_dma_submit 138.2/f 0.610/0.916/2.615 us 0.16%
```

- Mean pack 7.236 us/call versus submit 0.916 us/call (7.9×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.43%, corresponding to 60.981 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.791×; Show end: 0.748×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `rm_shader_draw` — 73.1% of the frame, 45.952 ms/f.
2. `rm_timeline_step` — 0.5% of the frame, 0.331 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.130 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Raymarch`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
