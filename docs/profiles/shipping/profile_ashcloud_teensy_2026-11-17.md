# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\ashcloud_ship.log`. Replaces `profile_ashcloud_teensy_2026-11-13.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh AshCloud profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:125896, data:163648, headers:8440; free:1733632 / RAM1: variables:314144, code:39096, padding:26440; free:144608 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 897–928 root cycles / 600 MHz = 1997622.540 us; measured wall sum 1997626 us, delta 1.73 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:39, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `f294f94a384f7bc31923bd25907792c91384443bcd7e6ca8bc8760c2aebd9b7a`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **85.379 ms**, spilled **1/960**. Setup frame 1 is **85.379 ms** and remains separate. Live show peak is **56.216 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 56.216, spilled 🟢 0/959 (0.00%).

README setup frame 1: 85.379 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.488 ms  37.493 Mcyc 100.0%
  pov_preserve_half      136.98 us   82.19 kcyc   0.2% x0.96875 141 us/call
  fx_shader_draw         40.758 ms  24.455 Mcyc  65.2%
    fx_prepare_column     1.758 ms   1.055 Mcyc   2.8% x144 12 us/call
  fx_palette_base_bake   578.12 us  346.87 kcyc   0.9%
    fx_palette_rgb_bake  451.30 us  270.78 kcyc   0.7% x1 451 us/call
  fx_prepare_frame        1.450 ms  869.78 kcyc   2.3% x1 1450 us/call
  fx_advance              1.104 ms  662.55 kcyc   1.8%
    palette_cycle_step    1.041 ms  624.64 kcyc   1.7% x1 1041 us/call
  fx_timeline_step       114.50 us   68.70 kcyc   0.2% x1 115 us/call
  canvas_clear            98.10 us   58.86 kcyc   0.2% x1 98 us/call
  canvas_buffer_wait     18.067 ms  10.840 Mcyc  28.9% x1 18067 us/call
```

Wall min/avg/max = 37.828/62.488/85.381 ms. `fx_shader_draw` averages 40.758 ms/f; total render averages 44.421 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 897–928)

```
frame                    62.426 ms  37.455 Mcyc 100.0%
  pov_preserve_half      136.18 us   81.71 kcyc   0.2% x1 136 us/call
  fx_shader_draw         48.116 ms  28.869 Mcyc  77.1%
    fx_prepare_column     1.738 ms   1.043 Mcyc   2.8% x144 12 us/call
  fx_palette_base_bake   938.16 us  562.89 kcyc   1.5%
    fx_palette_rgb_bake  798.56 us  479.14 kcyc   1.3% x1 799 us/call
  fx_prepare_frame        1.597 ms  958.29 kcyc   2.6% x1 1597 us/call
  fx_advance              1.108 ms  664.69 kcyc   1.8%
    palette_cycle_step    1.035 ms  620.77 kcyc   1.7% x1 1035 us/call
  fx_timeline_step       117.43 us   70.46 kcyc   0.2% x1 117 us/call
  canvas_clear            94.86 us   56.91 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     10.150 ms   6.090 Mcyc  16.3% x1 10150 us/call
```

Wall min/avg/max = 56.762/62.425/68.102 ms. `fx_shader_draw` averages 48.116 ms/f; total render averages 52.276 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.232 ms  37.339 Mcyc 100.0%
  pov_preserve_half      136.66 us   82.00 kcyc   0.2% x1 137 us/call
  fx_shader_draw         46.682 ms  28.009 Mcyc  75.0%
    fx_prepare_column     1.753 ms   1.052 Mcyc   2.8% x144 12 us/call
  fx_palette_base_bake   903.28 us  541.97 kcyc   1.5%
    fx_palette_rgb_bake  773.90 us  464.34 kcyc   1.2% x1 774 us/call
  fx_prepare_frame        1.592 ms  955.42 kcyc   2.6% x1 1592 us/call
  fx_advance              1.105 ms  663.09 kcyc   1.8%
    palette_cycle_step    1.036 ms  621.48 kcyc   1.7% x1 1036 us/call
  fx_timeline_step       115.82 us   69.49 kcyc   0.2% x1 116 us/call
  canvas_clear            94.26 us   56.56 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     11.435 ms   6.861 Mcyc  18.4% x1 11435 us/call
```

Wall min/avg/max = 57.304/62.232/67.052 ms. `fx_shader_draw` averages 46.682 ms/f; total render averages 50.797 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 897–928; CPU shares use its separately snapshotted ISR interval of 2001993 us.

```
isr_wake       1153.2/f 0.286/1.761/22.461 us 3.24%
  isr_pack       144.1/f 8.008/8.701/11.488 us 2.00%
  isr_dma_submit 144.1/f 0.656/0.926/8.971 us 0.21%
```

- Mean pack 8.701 us/call versus submit 0.926 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.24%, corresponding to 60.475 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.711×; Peak live show render: 0.836×; Show end: 0.813×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 77.1% of the frame, 48.116 ms/f.
2. `fx_prepare_column` — 2.8% of the frame, 1.738 ms/f.
3. `fx_prepare_frame` — 2.6% of the frame, 1.597 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
