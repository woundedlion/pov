# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\gridspace_ship.log`. Replaces `profile_gridspace_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | GridSpace 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh GridSpace profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:104616, data:152600, headers:9024; free:1765376 / RAM1: variables:314144, code:28456, padding:4312; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 353–384 root cycles / 600 MHz = 1996231.303 us; measured wall sum 1996227 us, delta 2.16 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:33, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `68aa9ea2eb300c0bf3e9f256c126327d20ed52a6281dc78bc2257ce92d2b2802`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **46.539 ms**, spilled **0/480**. Setup frame 1 is **46.539 ms** and remains separate. Live show peak is **32.794 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1096 frame records; the 616 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 32.794, spilled 🟢 0/479 (0.00%).

README setup frame 1: 46.539 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    60.696 ms  36.417 Mcyc 100.0%
  pov_preserve_half      138.77 us   83.26 kcyc   0.2% x0.96875 143 us/call
  fx_shader_draw         21.014 ms  12.609 Mcyc  34.6%
    fx_prepare_column    598.37 us  359.02 kcyc   1.0% x144 4 us/call
  fx_palette_base_bake   605.80 us  363.48 kcyc   1.0%
    fx_palette_rgb_bake  488.17 us  292.90 kcyc   0.8% x1 488 us/call
  fx_prepare_frame        3.523 ms   2.114 Mcyc   5.8% x1 3523 us/call
  fx_advance              1.031 ms  618.33 kcyc   1.7%
    palette_cycle_step   980.64 us  588.38 kcyc   1.6% x1 981 us/call
  fx_timeline_step       108.51 us   65.11 kcyc   0.2% x1 109 us/call
  canvas_clear            97.34 us   58.40 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     34.023 ms  20.414 Mcyc  56.1% x1 34023 us/call
```

Wall min/avg/max = 24.590/60.695/66.127 ms. `fx_shader_draw` averages 21.014 ms/f; total render averages 26.673 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 353–384)

```
frame                    62.382 ms  37.429 Mcyc 100.0%
  pov_preserve_half      144.02 us   86.41 kcyc   0.2% x1 144 us/call
  fx_shader_draw         21.130 ms  12.678 Mcyc  33.9%
    fx_prepare_column    603.87 us  362.32 kcyc   1.0% x144 4 us/call
  fx_palette_base_bake    1.294 ms  776.62 kcyc   2.1%
    fx_palette_rgb_bake   1.172 ms  703.42 kcyc   1.9% x1 1172 us/call
  fx_prepare_frame        3.656 ms   2.193 Mcyc   5.9% x1 3656 us/call
  fx_advance             995.37 us  597.22 kcyc   1.6%
    palette_cycle_step   946.35 us  567.81 kcyc   1.5% x1 946 us/call
  fx_timeline_step       115.67 us   69.41 kcyc   0.2% x1 116 us/call
  canvas_clear            92.48 us   55.49 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     34.803 ms  20.882 Mcyc  55.8% x1 34803 us/call
```

Wall min/avg/max = 55.271/62.382/69.506 ms. `fx_shader_draw` averages 21.130 ms/f; total render averages 27.579 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.378 ms  37.427 Mcyc 100.0%
  pov_preserve_half      145.09 us   87.06 kcyc   0.2% x1 145 us/call
  fx_shader_draw         20.753 ms  12.452 Mcyc  33.3%
    fx_prepare_column    601.21 us  360.73 kcyc   1.0% x144 4 us/call
  fx_palette_base_bake    1.170 ms  701.87 kcyc   1.9%
    fx_palette_rgb_bake   1.046 ms  627.42 kcyc   1.7% x1 1046 us/call
  fx_prepare_frame        3.633 ms   2.180 Mcyc   5.8% x1 3633 us/call
  fx_advance              1.026 ms  615.56 kcyc   1.6%
    palette_cycle_step   972.94 us  583.76 kcyc   1.6% x1 973 us/call
  fx_timeline_step       111.07 us   66.64 kcyc   0.2% x1 111 us/call
  canvas_clear            93.03 us   55.82 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     35.295 ms  21.177 Mcyc  56.6% x1 35295 us/call
```

Wall min/avg/max = 55.351/62.378/69.379 ms. `fx_shader_draw` averages 20.753 ms/f; total render averages 27.083 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 353–384; CPU shares use its separately snapshotted ISR interval of 2000269 us.

```
isr_wake       1152.2/f 0.288/1.666/21.775 us 3.07%
  isr_pack       144.0/f 7.423/8.291/10.868 us 1.91%
  isr_dma_submit 144.0/f 0.628/0.920/7.325 us 0.21%
```

- Mean pack 8.291 us/call versus submit 0.920 us/call (9.0×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.07%, corresponding to 60.581 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.427×; Peak live show render: 0.441×; Show end: 0.433×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 33.9% of the frame, 21.130 ms/f.
2. `fx_prepare_frame` — 5.9% of the frame, 3.656 ms/f.
3. `fx_palette_base_bake` — 2.1% of the frame, 1.294 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GridSpace`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
