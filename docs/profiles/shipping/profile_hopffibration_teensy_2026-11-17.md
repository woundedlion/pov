# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\hopffibration_ship.log`. Replaces `profile_hopffibration_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared scan-shapes, spherical-math and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HopfFibration profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:127680, data:145672, headers:8248; free:1750016 / RAM1: variables:312992, code:29160, padding:3608; free:178528 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 801–832 root cycles / 600 MHz = 2005423.490 us; measured wall sum 2005420 us, delta 1.74 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:45, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `32a4d29bfea69e1e9a2a3984ec26c11afd4f6f95b55b478154413801090c87c9`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **55.635 ms**, spilled **0/960**. Setup frame 1 is **0.672 ms** and remains separate. Live show peak is **55.635 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 55.635, spilled 🟢 0/959 (0.00%).

README setup frame 1: 0.672 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    59.294 ms  35.576 Mcyc 100.0%
  pov_preserve_half      139.59 us   83.76 kcyc   0.2% x0.96875 144 us/call
  hf_render_trails       10.762 ms   6.457 Mcyc  18.2% x1 10762 us/call
  hf_project_record      158.43 us   95.06 kcyc   0.3% x1 158 us/call
  hf_advance_tumble       10.29 us    6.18 kcyc   0.0% x1 10 us/call
  hf_timeline_step        62.48 us   37.49 kcyc   0.1% x1 62 us/call
  canvas_clear            95.49 us   57.29 kcyc   0.2% x1 95 us/call
  canvas_buffer_wait     48.045 ms  28.827 Mcyc  81.0% x1 48045 us/call
```

Wall min/avg/max = 0.672/59.293/69.575 ms. `hf_render_trails` averages 10.762 ms/f; total render averages 11.249 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 801–832)

```
frame                    62.669 ms  37.602 Mcyc 100.0%
  pov_preserve_half      136.28 us   81.77 kcyc   0.2% x1 136 us/call
  hf_render_trails       39.601 ms  23.761 Mcyc  63.2% x1 39601 us/call
  hf_project_record      187.21 us  112.33 kcyc   0.3% x1 187 us/call
  hf_advance_tumble       11.88 us    7.13 kcyc   0.0% x1 12 us/call
  hf_timeline_step        80.65 us   48.39 kcyc   0.1% x1 81 us/call
  canvas_clear            93.96 us   56.38 kcyc   0.1% x1 94 us/call
  canvas_buffer_wait     22.537 ms  13.522 Mcyc  36.0% x1 22537 us/call
```

Wall min/avg/max = 32.000/62.669/94.770 ms. `hf_render_trails` averages 39.601 ms/f; total render averages 40.132 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.836 ms  37.702 Mcyc 100.0%
  pov_preserve_half      139.76 us   83.85 kcyc   0.2% x1 140 us/call
  hf_render_trails       28.632 ms  17.179 Mcyc  45.6% x1 28632 us/call
  hf_project_record      187.24 us  112.34 kcyc   0.3% x1 187 us/call
  hf_advance_tumble        9.84 us    5.90 kcyc   0.0% x1 10 us/call
  hf_timeline_step        70.00 us   42.00 kcyc   0.1% x1 70 us/call
  canvas_clear            92.69 us   55.61 kcyc   0.1% x1 93 us/call
  canvas_buffer_wait     33.685 ms  20.211 Mcyc  53.6% x1 33685 us/call
```

Wall min/avg/max = 37.435/62.836/88.323 ms. `hf_render_trails` averages 28.632 ms/f; total render averages 29.151 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 801–832; CPU shares use its separately snapshotted ISR interval of 2008964 us.

```
isr_wake       1157.2/f 0.268/1.566/29.255 us 2.88%
  isr_pack       144.7/f 5.591/7.375/11.133 us 1.69%
  isr_dma_submit 144.7/f 0.575/0.913/10.698 us 0.21%
```

- Mean pack 7.375 us/call versus submit 0.913 us/call (8.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.88%, corresponding to 60.700 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.180×; Peak live show render: 0.642×; Show end: 0.466×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `hf_render_trails` — 63.2% of the frame, 39.601 ms/f.
2. `hf_project_record` — 0.3% of the frame, 0.187 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.136 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared scan-shapes, spherical-math and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HopfFibration`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
