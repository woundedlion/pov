# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\petalflow_ship.log`. Replaces `profile_petalflow_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared scan-shapes, spherical-math and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh PetalFlow profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:111400, data:145760, headers:9080; free:1765376 / RAM1: variables:312992, code:29688, padding:3080; free:178528 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 289–320 root cycles / 600 MHz = 2001151.367 us; measured wall sum 2001153 us, delta 0.82 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:50, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `c003f11aebc8b674ad47e13a86cc11b72de3341969e28cf5085b3068bab01f71`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **25.420 ms**, spilled **0/480**. Setup frame 1 is **19.712 ms** and remains separate. Live show peak is **25.420 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 1097 frame records; the 617 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 25.420, spilled 🟢 0/479 (0.00%).

README setup frame 1: 19.712 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    59.690 ms  35.814 Mcyc 100.0%
  pov_preserve_half      132.58 us   79.55 kcyc   0.2% x0.96875 137 us/call
  pf_draw_rings          18.855 ms  11.313 Mcyc  31.6%
    pf_ring_scan         17.855 ms  10.713 Mcyc  29.9% x23.3438 765 us/call
    pf_ring_build        854.19 us  512.52 kcyc   1.4% x23.3438 37 us/call
  pf_timeline_step        53.26 us   31.95 kcyc   0.1% x1 53 us/call
  canvas_clear            89.65 us   53.79 kcyc   0.2% x1 90 us/call
  canvas_buffer_wait     40.536 ms  24.321 Mcyc  67.9% x1 40536 us/call
```

Wall min/avg/max = 17.676/59.689/67.138 ms. `pf_draw_rings` averages 18.855 ms/f; total render averages 19.154 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 289–320)

```
frame                    62.536 ms  37.522 Mcyc 100.0%
  pov_preserve_half      134.87 us   80.92 kcyc   0.2% x1 135 us/call
  pf_draw_rings          23.108 ms  13.865 Mcyc  37.0%
    pf_ring_scan         22.099 ms  13.260 Mcyc  35.3% x23.3438 947 us/call
    pf_ring_build        854.48 us  512.69 kcyc   1.4% x23.3438 37 us/call
  pf_timeline_step        52.14 us   31.28 kcyc   0.1% x1 52 us/call
  canvas_clear            84.95 us   50.97 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     39.132 ms  23.479 Mcyc  62.6% x1 39132 us/call
```

Wall min/avg/max = 58.783/62.536/65.740 ms. `pf_draw_rings` averages 23.108 ms/f; total render averages 23.404 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.401 ms  37.441 Mcyc 100.0%
  pov_preserve_half      139.76 us   83.85 kcyc   0.2% x1 140 us/call
  pf_draw_rings          19.849 ms  11.909 Mcyc  31.8%
    pf_ring_scan         18.851 ms  11.311 Mcyc  30.2% x23.3438 808 us/call
    pf_ring_build        850.03 us  510.02 kcyc   1.4% x23.3438 36 us/call
  pf_timeline_step        53.20 us   31.92 kcyc   0.1% x1 53 us/call
  canvas_clear            84.76 us   50.86 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     42.251 ms  25.351 Mcyc  67.7% x1 42251 us/call
```

Wall min/avg/max = 60.421/62.401/64.890 ms. `pf_draw_rings` averages 19.849 ms/f; total render averages 20.150 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 289–320; CPU shares use its separately snapshotted ISR interval of 2003490 us.

```
isr_wake       1154.1/f 0.271/1.441/18.388 us 2.65%
  isr_pack       144.2/f 5.596/6.440/9.240 us 1.48%
  isr_dma_submit 144.2/f 0.580/0.908/8.023 us 0.20%
```

- Mean pack 6.440 us/call versus submit 0.908 us/call (7.1×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.65%, corresponding to 60.844 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.306×; Peak live show render: 0.374×; Show end: 0.322×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `pf_draw_rings` — 37.0% of the frame, 23.108 ms/f.
2. `pf_ring_scan` — 35.3% of the frame, 22.099 ms/f.
3. `pf_ring_build` — 1.4% of the frame, 0.854 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared scan-shapes, spherical-math and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=PetalFlow`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
