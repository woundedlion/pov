# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\bzreactiondiffusion_ship.log`. Replaces `profile_bzreactiondiffusion_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the effect-local physics/shading functions and shared scan/filter kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 130 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh BZReactionDiffusion profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:68512, data:156524, headers:8436; free:1798144 / RAM1: variables:314048, code:16296, padding:16472; free:177472 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 577–608 root cycles / 600 MHz = 1998386.767 us; measured wall sum 1998389 us, delta 1.12 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:32, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `19dd6675e85d2aea03ecbbe5d13c69965ebd22e0f7b26a486095207f3eda1767`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **76.068 ms**, spilled **1/960**. Setup frame 1 is **76.068 ms** and remains separate. Live show peak is **44.329 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 2056 frame records; the 1096 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 44.329, spilled 🟢 0/959 (0.00%).

README setup frame 1: 76.068 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    62.304 ms  37.382 Mcyc 100.0%
  pov_preserve_half      133.20 us   79.92 kcyc   0.2% x0.96875 138 us/call
  bz_render              44.024 ms  26.414 Mcyc  70.7%
    bz_raster            39.766 ms  23.860 Mcyc  63.8% x2 19883 us/call
    bz_physics            4.164 ms   2.499 Mcyc   6.7% x2 2082 us/call
    rd_lattice_basis      90.02 us   54.01 kcyc   0.1% x1 90 us/call
  rd_timeline_step        32.58 us   19.55 kcyc   0.1% x1 33 us/call
  canvas_clear            89.70 us   53.82 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     18.010 ms  10.806 Mcyc  28.9% x1 18010 us/call
```

Wall min/avg/max = 42.051/62.303/76.069 ms. `bz_render` averages 44.024 ms/f; total render averages 44.294 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 577–608)

```
frame                    62.450 ms  37.470 Mcyc 100.0%
  pov_preserve_half      138.17 us   82.90 kcyc   0.2% x1 138 us/call
  bz_render              43.718 ms  26.231 Mcyc  70.0%
    bz_raster            39.499 ms  23.700 Mcyc  63.2% x2 19750 us/call
    bz_physics            4.127 ms   2.476 Mcyc   6.6% x2 2063 us/call
    rd_lattice_basis      89.43 us   53.66 kcyc   0.1% x1 89 us/call
  rd_timeline_step        28.41 us   17.05 kcyc   0.0% x1 28 us/call
  canvas_clear            85.79 us   51.48 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     18.465 ms  11.079 Mcyc  29.6% x1 18465 us/call
```

Wall min/avg/max = 61.638/62.449/63.106 ms. `bz_render` averages 43.718 ms/f; total render averages 43.984 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.454 ms  37.472 Mcyc 100.0%
  pov_preserve_half      134.50 us   80.70 kcyc   0.2% x1 135 us/call
  bz_render              43.589 ms  26.154 Mcyc  69.8%
    bz_raster            39.354 ms  23.612 Mcyc  63.0% x2 19677 us/call
    bz_physics            4.143 ms   2.486 Mcyc   6.6% x2 2072 us/call
    rd_lattice_basis      89.43 us   53.66 kcyc   0.1% x1 89 us/call
  rd_timeline_step        28.07 us   16.84 kcyc   0.0% x1 28 us/call
  canvas_clear            85.78 us   51.47 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     18.602 ms  11.161 Mcyc  29.8% x1 18602 us/call
```

Wall min/avg/max = 61.781/62.454/63.012 ms. `bz_render` averages 43.589 ms/f; total render averages 43.852 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 577–608; CPU shares use its separately snapshotted ISR interval of 2000015 us.

```
isr_wake       1152.1/f 0.268/1.733/13.813 us 3.19%
  isr_pack       144.0/f 8.125/8.701/11.743 us 2.00%
  isr_dma_submit 144.0/f 0.588/0.913/1.201 us 0.21%
```

- Mean pack 8.701 us/call versus submit 0.913 us/call (9.5×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.19%, corresponding to 60.506 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.709×; Peak live show render: 0.704×; Show end: 0.702×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `bz_render` — 70.0% of the frame, 43.718 ms/f.
2. `bz_raster` — 63.2% of the frame, 39.499 ms/f.
3. `bz_physics` — 6.6% of the frame, 4.127 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the effect-local physics/shading functions and shared scan/filter kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=BZReactionDiffusion`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
