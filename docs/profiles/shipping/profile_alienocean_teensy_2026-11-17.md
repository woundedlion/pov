# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\alienocean_ship.log`. Replaces `profile_alienocean_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | AlienOcean 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh AlienOcean profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:104688, data:152784, headers:8768; free:1765376 / RAM1: variables:314144, code:29128, padding:3640; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 385–416 root cycles / 600 MHz = 1997077.862 us; measured wall sum 1997077 us, delta 0.43 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:31, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `2877aad483a9eba596db8d9c14991dbc74786d0b953c4fd5d0e7662129034b57`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **49.224 ms**, spilled **0/480**. Setup frame 1 is **49.224 ms** and remains separate. Live show peak is **27.306 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1097 frame records; the 617 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 27.306, spilled 🟢 0/479 (0.00%).

README setup frame 1: 49.224 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    60.779 ms  36.468 Mcyc 100.0%
  pov_preserve_half      139.01 us   83.40 kcyc   0.2% x0.96875 143 us/call
  fx_shader_draw         22.854 ms  13.712 Mcyc  37.6%
    fx_prepare_column    861.25 us  516.75 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake   560.89 us  336.53 kcyc   0.9%
    fx_palette_rgb_bake  443.27 us  265.96 kcyc   0.7% x1 443 us/call
  fx_prepare_frame       816.38 us  489.83 kcyc   1.3% x1 816 us/call
  fx_advance             963.88 us  578.33 kcyc   1.6%
    palette_cycle_step   924.14 us  554.48 kcyc   1.5% x1 924 us/call
  fx_timeline_step        78.70 us   47.22 kcyc   0.1% x1 79 us/call
  canvas_clear            96.97 us   58.18 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     35.129 ms  21.077 Mcyc  57.8% x1 35129 us/call
```

Wall min/avg/max = 23.213/60.779/66.209 ms. `fx_shader_draw` averages 22.854 ms/f; total render averages 25.651 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 385–416)

```
frame                    62.409 ms  37.445 Mcyc 100.0%
  pov_preserve_half      140.56 us   84.34 kcyc   0.2% x1 141 us/call
  fx_shader_draw         22.173 ms  13.304 Mcyc  35.5%
    fx_prepare_column    853.37 us  512.02 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake   915.07 us  549.04 kcyc   1.5%
    fx_palette_rgb_bake  785.83 us  471.50 kcyc   1.3% x1 786 us/call
  fx_prepare_frame       760.42 us  456.25 kcyc   1.2% x1 760 us/call
  fx_advance              1.060 ms  635.88 kcyc   1.7%
    palette_cycle_step    1.023 ms  613.68 kcyc   1.6% x1 1023 us/call
  fx_timeline_step        79.21 us   47.52 kcyc   0.1% x1 79 us/call
  canvas_clear            93.00 us   55.80 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     37.042 ms  22.225 Mcyc  59.4% x1 37042 us/call
```

Wall min/avg/max = 58.736/62.408/66.140 ms. `fx_shader_draw` averages 22.173 ms/f; total render averages 25.367 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.391 ms  37.434 Mcyc 100.0%
  pov_preserve_half      140.68 us   84.41 kcyc   0.2% x1 141 us/call
  fx_shader_draw         22.189 ms  13.314 Mcyc  35.6%
    fx_prepare_column    854.24 us  512.55 kcyc   1.4% x144 6 us/call
  fx_palette_base_bake   796.60 us  477.96 kcyc   1.3%
    fx_palette_rgb_bake  667.76 us  400.66 kcyc   1.1% x1 668 us/call
  fx_prepare_frame       756.26 us  453.76 kcyc   1.2% x1 756 us/call
  fx_advance              1.057 ms  634.43 kcyc   1.7%
    palette_cycle_step    1.022 ms  613.31 kcyc   1.6% x1 1022 us/call
  fx_timeline_step        76.16 us   45.70 kcyc   0.1% x1 76 us/call
  canvas_clear            92.85 us   55.71 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     37.132 ms  22.279 Mcyc  59.5% x1 37132 us/call
```

Wall min/avg/max = 58.668/62.390/66.119 ms. `fx_shader_draw` averages 22.189 ms/f; total render averages 25.259 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 385–416; CPU shares use its separately snapshotted ISR interval of 2000044 us.

```
isr_wake       1152.1/f 0.285/1.720/14.448 us 3.17%
  isr_pack       144.0/f 6.320/8.608/11.198 us 1.98%
  isr_dma_submit 144.0/f 0.650/0.915/4.398 us 0.21%
```

- Mean pack 8.608 us/call versus submit 0.915 us/call (9.4×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.17%, corresponding to 60.519 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.410×; Peak live show render: 0.406×; Show end: 0.404×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 35.5% of the frame, 22.173 ms/f.
2. `fx_advance` — 1.7% of the frame, 1.060 ms/f.
3. `palette_cycle_step` — 1.6% of the frame, 1.023 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienOcean`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
