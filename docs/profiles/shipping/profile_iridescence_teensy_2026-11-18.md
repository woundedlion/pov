# Iridescence on-device profile — Teensy 4.0, segmented mode (2026-11-18, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\river59-evidence\review-finishiri.log`. Replaces `profile_iridescence_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/shader.h` and the pullback runtime/palette kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Iridescence 288×144, single-entry playlist, source `62520ec25e068949eeb144e34cfd023ef05f4daf` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 34 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Iridescence profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:110776, data:157204, headers:8500; free:1755136 / RAM1: variables:314144, code:30344, padding:2424; free:177376 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1089024, data:429272, headers:8480; free:504840 / RAM1: variables:313760, code:190888, padding:5720; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 257–288 root cycles / 600 MHz = 1996218.633 us; measured wall sum 1996221 us, delta 1.19 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-18 18:14, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `dfae460a280d445d7fb834a6b5268e92439176fbfa906b7ef99726b6a324979d`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Shared phase shading

The pre-review build peaked at 45.515 ms on COM3. The reviewed build peaked at 44.915 ms over 479 live frames with 0 display deadline misses. Its 480 native frames matched the earlier RGB16 capture byte for byte.

### Earlier placement measurements

The previous shipping renderer peaked at 47.242 ms on COM3. Sharing more code across the show brought that peak to 48.252 ms. The phase stack now executes its pixel shader and signed-byte cube blending in Teensy instruction RAM. A separate run of this placement peaked at 45.513 ms, and this qualification peaked at 45.515 ms with zero missed live frames. The same 480 native quadrants matched the earlier RGB16 capture byte for byte.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **83.018 ms**, spilled **1/480**. Setup frame 1 is **83.018 ms** and remains separate. Live show peak is **44.915 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 520 frame records; the 40 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 44.915, spilled 🟢 0/479 (0.00%).

README setup frame 1: 83.018 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.484 ms  37.491 Mcyc 100.0%
  pov_preserve_half      132.93 us   79.76 kcyc   0.2% x0.96875 137 us/call
  fx_shader_draw         42.644 ms  25.587 Mcyc  68.2%
    fx_prepare_column    984.26 us  590.56 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake    93.20 us   55.92 kcyc   0.1% x1 93 us/call
  fx_prepare_frame        2.526 ms   1.516 Mcyc   4.0% x1 2526 us/call
  fx_advance              24.64 us   14.78 kcyc   0.0% x1 25 us/call
  fx_timeline_step        12.05 us    7.23 kcyc   0.0% x1 12 us/call
  canvas_clear            91.02 us   54.61 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     16.744 ms  10.046 Mcyc  26.8% x1 16744 us/call
```

Wall min/avg/max = 44.416/62.484/83.018 ms. `fx_shader_draw` averages 42.644 ms/f; total render averages 45.740 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 257–288)

```
frame                    62.382 ms  37.429 Mcyc 100.0%
  pov_preserve_half      137.84 us   82.71 kcyc   0.2% x1 138 us/call
  fx_shader_draw         41.471 ms  24.882 Mcyc  66.5%
    fx_prepare_column     1.001 ms  600.41 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake    91.50 us   54.90 kcyc   0.1% x1 92 us/call
  fx_prepare_frame        2.455 ms   1.473 Mcyc   3.9% x1 2455 us/call
  fx_advance              24.40 us   14.64 kcyc   0.0% x1 24 us/call
  fx_timeline_step        10.30 us    6.18 kcyc   0.0% x1 10 us/call
  canvas_clear            86.69 us   52.02 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     17.891 ms  10.735 Mcyc  28.7% x1 17891 us/call
```

Wall min/avg/max = 61.575/62.381/63.102 ms. `fx_shader_draw` averages 41.471 ms/f; total render averages 44.491 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.391 ms  37.435 Mcyc 100.0%
  pov_preserve_half      140.52 us   84.31 kcyc   0.2% x1 141 us/call
  fx_shader_draw         41.357 ms  24.814 Mcyc  66.3%
    fx_prepare_column     1.000 ms  600.07 kcyc   1.6% x144 7 us/call
  fx_palette_base_bake    92.87 us   55.72 kcyc   0.1% x1 93 us/call
  fx_prepare_frame        2.449 ms   1.469 Mcyc   3.9% x1 2449 us/call
  fx_advance              24.14 us   14.48 kcyc   0.0% x1 24 us/call
  fx_timeline_step        12.07 us    7.24 kcyc   0.0% x1 12 us/call
  canvas_clear            86.74 us   52.04 kcyc   0.1% x1 87 us/call
  canvas_buffer_wait     18.015 ms  10.809 Mcyc  28.9% x1 18015 us/call
```

Wall min/avg/max = 62.055/62.391/62.766 ms. `fx_shader_draw` averages 41.357 ms/f; total render averages 44.376 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 257–288; CPU shares use its separately snapshotted ISR interval of 1999986 us.

```
isr_wake       1152.1/f 0.473/1.743/20.886 us 3.21%
  isr_pack       144.0/f 6.795/7.516/10.346 us 1.73%
  isr_dma_submit 144.0/f 0.633/0.915/2.451 us 0.21%
```

- Mean pack 7.516 us/call versus submit 0.915 us/call (8.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.21%, corresponding to 60.494 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.732×; Peak live show render: 0.712×; Show end: 0.710×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fx_shader_draw` — 66.5% of the frame, 41.471 ms/f.
2. `fx_prepare_frame` — 3.9% of the frame, 2.455 ms/f.
3. `fx_prepare_column` — 1.6% of the frame, 1.001 ms/f.

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
