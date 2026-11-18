# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\alienbrain_ship.log`. Replaces `profile_alienbrain_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM4, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | AlienBrain 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 240 s / 3840 nominal frames at 16 fps; window 16. Raw capture 285 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 bash tools/profile_one.sh AlienBrain profile 243 16`; analyze only the 240 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:107104, data:153192, headers:9016; free:1762304 / RAM1: variables:314144, code:28888, padding:3880; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 2033–2048 root cycles / 600 MHz = 997774.978 us; measured wall sum 997776 us, delta 1.02 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:13, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `9f72e533e1423c1495dd49f4b4c334437ec3aa9f7f6887e868919525d4556d86`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 240 s (3840 nominal frames at 16 fps). Show-frame peak render is **57.420 ms**, spilled **0/3840**. Setup frame 1 is **57.420 ms** and remains separate. Live show peak is **37.447 ms**, spilled **0/3839**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–3840. The raw log retains 4536 frame records; the 696 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 37.447 (4), spilled 🟢 0/3839 (0.00%).

README setup frame 1: 57.420 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    59.916 ms  35.950 Mcyc 100.0%
  pov_preserve_half      132.50 us   79.50 kcyc   0.2% x0.9375 141 us/call
  fx_shader_draw         27.915 ms  16.749 Mcyc  46.6%
    fx_prepare_column    765.66 us  459.39 kcyc   1.3% x144 5 us/call
  fx_palette_base_bake   634.21 us  380.53 kcyc   1.1%
    fx_palette_rgb_bake  502.34 us  301.41 kcyc   0.8% x1 502 us/call
  fx_prepare_frame       794.37 us  476.62 kcyc   1.3% x1 794 us/call
  fx_advance             956.81 us  574.08 kcyc   1.6%
    palette_cycle_step   886.73 us  532.04 kcyc   1.5% x1 887 us/call
  fx_timeline_step        99.57 us   59.74 kcyc   0.2% x1 100 us/call
  canvas_clear            99.66 us   59.79 kcyc   0.2% x1 100 us/call
  canvas_buffer_wait     29.130 ms  17.478 Mcyc  48.6% x1 29130 us/call
```

Wall min/avg/max = 27.933/59.916/65.647 ms. `fx_shader_draw` averages 27.915 ms/f; total render averages 30.786 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 2033–2048)

```
frame                    62.361 ms  37.417 Mcyc 100.0%
  pov_preserve_half      141.10 us   84.66 kcyc   0.2% x1 141 us/call
  fx_shader_draw         28.526 ms  17.116 Mcyc  45.7%
    fx_prepare_column     1.786 ms   1.072 Mcyc   2.9% x144 12 us/call
  fx_palette_base_bake    1.648 ms  988.79 kcyc   2.6%
    fx_palette_rgb_bake   1.513 ms  908.05 kcyc   2.4% x1 1513 us/call
  fx_prepare_frame       761.35 us  456.81 kcyc   1.2% x1 761 us/call
  fx_advance             977.91 us  586.74 kcyc   1.6%
    palette_cycle_step   916.10 us  549.66 kcyc   1.5% x1 916 us/call
  fx_timeline_step       181.70 us  109.02 kcyc   0.3% x1 182 us/call
  canvas_clear            90.94 us   54.57 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     29.887 ms  17.932 Mcyc  47.9% x1 29887 us/call
```

Wall min/avg/max = 54.964/62.361/69.759 ms. `fx_shader_draw` averages 28.526 ms/f; total render averages 32.474 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 3521–3536)

```
frame                    62.361 ms  37.416 Mcyc 100.0%
  pov_preserve_half      141.03 us   84.62 kcyc   0.2% x1 141 us/call
  fx_shader_draw         27.504 ms  16.502 Mcyc  44.1%
    fx_prepare_column    776.68 us  466.01 kcyc   1.2% x144 5 us/call
  fx_palette_base_bake    1.040 ms  624.22 kcyc   1.7%
    fx_palette_rgb_bake  905.68 us  543.41 kcyc   1.5% x1 906 us/call
  fx_prepare_frame       761.66 us  456.99 kcyc   1.2% x1 762 us/call
  fx_advance             964.79 us  578.87 kcyc   1.5%
    palette_cycle_step   892.58 us  535.55 kcyc   1.4% x1 893 us/call
  fx_timeline_step        99.37 us   59.62 kcyc   0.2% x1 99 us/call
  canvas_clear            91.41 us   54.84 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     31.606 ms  18.964 Mcyc  50.7% x1 31606 us/call
```

Wall min/avg/max = 58.738/62.360/65.926 ms. `fx_shader_draw` averages 27.504 ms/f; total render averages 30.754 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 897–912)

```
frame                    62.400 ms  37.440 Mcyc 100.0%
  pov_preserve_half      141.19 us   84.71 kcyc   0.2% x1 141 us/call
  fx_shader_draw         28.777 ms  17.266 Mcyc  46.1%
    fx_prepare_column     1.793 ms   1.076 Mcyc   2.9% x144 12 us/call
  fx_palette_base_bake    2.120 ms   1.272 Mcyc   3.4%
    fx_palette_rgb_bake   1.983 ms   1.190 Mcyc   3.2% x1 1983 us/call
  fx_prepare_frame       659.04 us  395.42 kcyc   1.1% x1 659 us/call
  fx_advance             973.54 us  584.12 kcyc   1.6%
    palette_cycle_step   912.43 us  547.46 kcyc   1.5% x1 912 us/call
  fx_timeline_step       179.31 us  107.59 kcyc   0.3% x1 179 us/call
  canvas_clear            90.27 us   54.16 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     29.310 ms  17.586 Mcyc  47.0% x1 29310 us/call
```

Wall min/avg/max = 55.282/62.400/69.314 ms. `fx_shader_draw` averages 28.777 ms/f; total render averages 33.090 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 3825–3840)

```
frame                    62.445 ms  37.467 Mcyc 100.0%
  pov_preserve_half      139.56 us   83.74 kcyc   0.2% x1 140 us/call
  fx_shader_draw         26.324 ms  15.795 Mcyc  42.2%
    fx_prepare_column    969.49 us  581.70 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake    1.979 ms   1.188 Mcyc   3.2%
    fx_palette_rgb_bake   1.845 ms   1.107 Mcyc   3.0% x1 1845 us/call
  fx_prepare_frame       722.72 us  433.63 kcyc   1.2% x1 723 us/call
  fx_advance             974.59 us  584.76 kcyc   1.6%
    palette_cycle_step   905.58 us  543.35 kcyc   1.5% x1 906 us/call
  fx_timeline_step       121.18 us   72.71 kcyc   0.2% x1 121 us/call
  canvas_clear            89.98 us   53.99 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     31.948 ms  19.169 Mcyc  51.2% x1 31948 us/call
```

Wall min/avg/max = 55.736/62.445/69.006 ms. `fx_shader_draw` averages 26.324 ms/f; total render averages 30.497 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held windows use fixed authored parameters while the scene keeps animating. Qualified windows exclude incoming parameter blends, fade-in and fade-out, and mixed-owner windows.

The show-scoped data covers 4 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 4 · Alien Brain 4 | 36 | 3521–3536 | held | — | — | 27.504 | 30.754 | 16 | 35.404 | 0/1079 |
| 1 · Alien Brain | 36 | 337–352 | held | — | — | 27.371 | 30.878 | 16 | 35.683 | 0/602 |
| 2 · Alien Brain 2 | 36 | 1441–1456 | held | — | — | 27.015 | 29.817 | 16 | 36.703 | 0/1079 |
| 3 · Alien Brain 3 | 37 | 2593–2608 | held | — | — | 26.025 | 30.236 | 16 | 37.447 | 0/1079 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 2033–2048; CPU shares use its separately snapshotted ISR interval of 999569 us.

```
isr_wake       1151.6/f 0.328/1.678/27.680 us 3.09%
  isr_pack       143.9/f 7.140/8.155/10.811 us 1.87%
  isr_dma_submit 143.9/f 0.661/0.913/0.980 us 0.21%
```

- Mean pack 8.155 us/call versus submit 0.913 us/call (8.9×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.09%, corresponding to 60.569 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.493×; Peak live show render: 0.520×; Transition / non-hold regime: 0.529×; Show end: 0.488×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 45.7% of the frame, 28.526 ms/f.
2. `fx_prepare_column` — 2.9% of the frame, 1.786 ms/f.
3. `fx_palette_base_bake` — 2.6% of the frame, 1.648 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/shader.h` and the pullback runtime/palette kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienBrain`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
