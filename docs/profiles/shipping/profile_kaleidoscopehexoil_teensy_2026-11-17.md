# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\kaleidoscopehexoil_ship.log`. Replaces `profile_kaleidoscopehexoil_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 105 s / 1680 nominal frames at 16 fps; window 16. Raw capture 150 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeHexOil profile 108 16`; analyze only the 105 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:101032, data:157012, headers:8196; free:1765376 / RAM1: variables:314144, code:26760, padding:6008; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 849–864 root cycles / 600 MHz = 997098.198 us; measured wall sum 997101 us, delta 2.81 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:26, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `8acdac244714a5758aed1a3428609f4e6d30085d454acdd3d714dcae3ad7f81c`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 105 s (1680 nominal frames at 16 fps). Show-frame peak render is **68.013 ms**, spilled **1/1680**. Setup frame 1 is **68.013 ms** and remains separate. Live show peak is **48.354 ms**, spilled **0/1679**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1680. The raw log retains 2376 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 48.354 (2), spilled 🟢 0/1679 (0.00%).

README setup frame 1: 68.013 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    60.992 ms  36.595 Mcyc 100.0%
  pov_preserve_half      130.68 us   78.41 kcyc   0.2% x0.9375 139 us/call
  fx_shader_draw         35.421 ms  21.253 Mcyc  58.1%
    fx_prepare_column     1.234 ms  740.40 kcyc   2.0% x144 9 us/call
  fx_palette_base_bake   640.48 us  384.29 kcyc   1.1%
    fx_palette_rgb_bake  509.74 us  305.84 kcyc   0.8% x1 510 us/call
  fx_prepare_frame       748.72 us  449.23 kcyc   1.2% x1 749 us/call
  fx_advance              1.130 ms  678.12 kcyc   1.9%
    palette_cycle_step    1.061 ms  636.49 kcyc   1.7% x1 1061 us/call
  fx_timeline_step        96.39 us   57.84 kcyc   0.2% x1 96 us/call
  canvas_clear           100.24 us   60.14 kcyc   0.2% x1 100 us/call
  canvas_buffer_wait     22.506 ms  13.503 Mcyc  36.9% x1 22506 us/call
```

Wall min/avg/max = 35.142/60.991/68.016 ms. `fx_shader_draw` averages 35.421 ms/f; total render averages 38.486 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 849–864)

```
frame                    62.319 ms  37.391 Mcyc 100.0%
  pov_preserve_half      139.50 us   83.70 kcyc   0.2% x1 140 us/call
  fx_shader_draw         39.028 ms  23.417 Mcyc  62.6%
    fx_prepare_column     2.616 ms   1.570 Mcyc   4.2% x144 18 us/call
  fx_palette_base_bake    2.049 ms   1.230 Mcyc   3.3%
    fx_palette_rgb_bake   1.909 ms   1.145 Mcyc   3.1% x1 1909 us/call
  fx_prepare_frame       834.75 us  500.85 kcyc   1.3% x1 835 us/call
  fx_advance              1.108 ms  664.83 kcyc   1.8%
    palette_cycle_step    1.046 ms  627.52 kcyc   1.7% x1 1046 us/call
  fx_timeline_step       184.94 us  110.96 kcyc   0.3% x1 185 us/call
  canvas_clear            91.01 us   54.60 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     18.681 ms  11.209 Mcyc  30.0% x1 18681 us/call
```

Wall min/avg/max = 53.687/62.318/70.912 ms. `fx_shader_draw` averages 39.028 ms/f; total render averages 43.638 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 353–368)

```
frame                    62.338 ms  37.403 Mcyc 100.0%
  pov_preserve_half      136.87 us   82.12 kcyc   0.2% x1 137 us/call
  fx_shader_draw         36.721 ms  22.033 Mcyc  58.9%
    fx_prepare_column     1.397 ms  838.40 kcyc   2.2% x144 10 us/call
  fx_palette_base_bake    1.336 ms  801.37 kcyc   2.1%
    fx_palette_rgb_bake   1.196 ms  717.61 kcyc   1.9% x1 1196 us/call
  fx_prepare_frame       949.71 us  569.83 kcyc   1.5% x1 950 us/call
  fx_advance              1.093 ms  655.97 kcyc   1.8%
    palette_cycle_step    1.024 ms  614.49 kcyc   1.6% x1 1024 us/call
  fx_timeline_step        95.73 us   57.44 kcyc   0.2% x1 96 us/call
  canvas_clear            91.35 us   54.81 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     21.704 ms  13.022 Mcyc  34.8% x1 21704 us/call
```

Wall min/avg/max = 58.477/62.337/66.153 ms. `fx_shader_draw` averages 36.721 ms/f; total render averages 40.634 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1665–1680)

```
frame                    62.417 ms  37.450 Mcyc 100.0%
  pov_preserve_half      138.37 us   83.02 kcyc   0.2% x1 138 us/call
  fx_shader_draw         33.717 ms  20.230 Mcyc  54.0%
    fx_prepare_column     1.311 ms  786.57 kcyc   2.1% x144 9 us/call
  fx_palette_base_bake   888.12 us  532.87 kcyc   1.4%
    fx_palette_rgb_bake  752.11 us  451.27 kcyc   1.2% x1 752 us/call
  fx_prepare_frame       806.31 us  483.78 kcyc   1.3% x1 806 us/call
  fx_advance              1.087 ms  652.07 kcyc   1.7%
    palette_cycle_step    1.015 ms  609.27 kcyc   1.6% x1 1015 us/call
  fx_timeline_step        99.54 us   59.73 kcyc   0.2% x1 100 us/call
  canvas_clear            91.24 us   54.74 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     25.386 ms  15.232 Mcyc  40.7% x1 25386 us/call
```

Wall min/avg/max = 60.904/62.417/63.732 ms. `fx_shader_draw` averages 33.717 ms/f; total render averages 37.031 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 2 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · Kaleidoscope Hex Oil | 36 | 353–368 | held | — | — | 36.721 | 40.634 | 16 | 46.472 | 0/600 |
| 2 · Kaleidoscope Hex Oil 2 | 36 | 1297–1312 | held | — | — | 35.664 | 38.709 | 16 | 48.354 | 0/1079 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 849–864; CPU shares use its separately snapshotted ISR interval of 999723 us.

```
isr_wake       1151.8/f 0.395/1.733/24.515 us 3.19%
  isr_pack       144.0/f 7.368/8.428/11.105 us 1.94%
  isr_dma_submit 144.0/f 0.576/0.921/5.708 us 0.21%
```

- Mean pack 8.428 us/call versus submit 0.921 us/call (9.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.19%, corresponding to 60.506 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.616×; Peak live show render: 0.698×; Costliest qualified preset: 0.650×; Show end: 0.592×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 62.6% of the frame, 39.028 ms/f.
2. `fx_prepare_column` — 4.2% of the frame, 2.616 ms/f.
3. `fx_palette_base_bake` — 3.3% of the frame, 2.049 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexOil`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
