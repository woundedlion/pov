# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\shared-mod-evidence\iri-paired-bands.log`. Replaces `profile_iridescence_teensy_2026-11-17.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Iridescence 288×144, single-entry playlist, source `e2e5353c810d11b962dbaa945932da6a9af06fce` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Iridescence profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:114816, data:157296, headers:8464; free:1751040 / RAM1: variables:314144, code:26792, padding:5976; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1172064, data:430244, headers:8436; free:420872 / RAM1: variables:313760, code:175528, padding:21080; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 65–96 root cycles / 600 MHz = 1997418.620 us; measured wall sum 1997420 us, delta 0.69 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 11:18, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `4e3e0172c8322c64c30825cd94383dd101eea73bd3195e8f071027098547b311`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Shared modulation comparison

On the same COM3 Teensy, the original renderer peaked at 48.106 ms. The shared renderer peaked at 47.242 ms. Both completed all 479 live frames within 62.5 ms. Setup fell from 91.760 to 88.278 ms. A 480-frame native comparison matched every output pixel (RGB16 SHA-256 `cd7cb04a2227559ed7f7b807d9028c491c415f6bc6cd7b76c8d5eb6b177a4be7`). The byte-cube evaluator prepares one set of tap coordinates per pixel for the ordered layers and reuses equal previous/current corner values.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **88.278 ms**, spilled **1/480**. Setup frame 1 is **88.278 ms** and remains separate. Live show peak is **47.242 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 520 frame records; the 40 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 47.242, spilled 🟢 0/479 (0.00%).

README setup frame 1: 88.278 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.683 ms  37.610 Mcyc 100.0%
  pov_preserve_half      133.51 us   80.10 kcyc   0.2% x0.96875 138 us/call
  fx_shader_draw         42.474 ms  25.485 Mcyc  67.8%
    fx_prepare_column     1.241 ms  744.65 kcyc   2.0% x144 9 us/call
  fx_palette_base_bake    75.84 us   45.51 kcyc   0.1% x1 76 us/call
  fx_prepare_frame        3.233 ms   1.940 Mcyc   5.2% x1 3233 us/call
  fx_advance             987.54 us  592.52 kcyc   1.6%
    palette_cycle_step   966.92 us  580.15 kcyc   1.5% x1 967 us/call
  fx_timeline_step        14.36 us    8.62 kcyc   0.0% x1 14 us/call
  canvas_clear            92.76 us   55.66 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     15.431 ms   9.259 Mcyc  24.6% x1 15431 us/call
```

Wall min/avg/max = 45.008/62.682/88.280 ms. `fx_shader_draw` averages 42.474 ms/f; total render averages 47.251 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 65–96)

```
frame                    62.419 ms  37.452 Mcyc 100.0%
  pov_preserve_half      137.62 us   82.57 kcyc   0.2% x1 138 us/call
  fx_shader_draw         41.274 ms  24.765 Mcyc  66.1%
    fx_prepare_column     1.248 ms  749.10 kcyc   2.0% x144 9 us/call
  fx_palette_base_bake    73.35 us   44.01 kcyc   0.1% x1 73 us/call
  fx_prepare_frame        3.059 ms   1.835 Mcyc   4.9% x1 3059 us/call
  fx_advance             943.35 us  566.01 kcyc   1.5%
    palette_cycle_step   923.28 us  553.97 kcyc   1.5% x1 923 us/call
  fx_timeline_step        15.02 us    9.01 kcyc   0.0% x1 15 us/call
  canvas_clear            87.59 us   52.55 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     16.592 ms   9.955 Mcyc  26.6% x1 16592 us/call
```

Wall min/avg/max = 60.078/62.419/64.733 ms. `fx_shader_draw` averages 41.274 ms/f; total render averages 45.827 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.422 ms  37.453 Mcyc 100.0%
  pov_preserve_half      136.56 us   81.94 kcyc   0.2% x1 137 us/call
  fx_shader_draw         41.198 ms  24.719 Mcyc  66.0%
    fx_prepare_column     1.249 ms  749.23 kcyc   2.0% x144 9 us/call
  fx_palette_base_bake    74.07 us   44.44 kcyc   0.1% x1 74 us/call
  fx_prepare_frame        3.060 ms   1.836 Mcyc   4.9% x1 3060 us/call
  fx_advance             813.28 us  487.97 kcyc   1.3%
    palette_cycle_step   793.19 us  475.91 kcyc   1.3% x1 793 us/call
  fx_timeline_step        15.12 us    9.07 kcyc   0.0% x1 15 us/call
  canvas_clear            88.33 us   53.00 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     16.800 ms  10.080 Mcyc  26.9% x1 16800 us/call
```

Wall min/avg/max = 60.677/62.421/64.034 ms. `fx_shader_draw` averages 41.198 ms/f; total render averages 45.621 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 65–96; CPU shares use its separately snapshotted ISR interval of 2000094 us.

```
isr_wake       1152.1/f 0.448/1.760/21.126 us 3.24%
  isr_pack       144.0/f 6.821/7.526/10.083 us 1.73%
  isr_dma_submit 144.0/f 0.571/0.916/9.725 us 0.21%
```

- Mean pack 7.526 us/call versus submit 0.916 us/call (8.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.24%, corresponding to 60.475 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.756×; Peak live show render: 0.733×; Show end: 0.730×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 66.1% of the frame, 41.274 ms/f.
2. `fx_prepare_frame` — 4.9% of the frame, 3.059 ms/f.
3. `fx_prepare_column` — 2.0% of the frame, 1.248 ms/f.

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
