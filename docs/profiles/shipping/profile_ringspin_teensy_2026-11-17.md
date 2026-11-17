# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\ringspin_ship.log`. Replaces `profile_ringspin_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared Plot raster, column-walk and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh RingSpin profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:86880, data:148020, headers:8812; free:1787904 / RAM1: variables:313952, code:27384, padding:5384; free:177568 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 737–768 root cycles / 600 MHz = 2003065.542 us; measured wall sum 2003064 us, delta 0.77 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:54, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `fc275bfc410133f6f3837f9b48e39338c21762a2c85590cf996200003e9368a4`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **51.117 ms**, spilled **0/960**. Setup frame 1 is **3.279 ms** and remains separate. Live show peak is **51.117 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 51.117, spilled 🟢 0/959 (0.00%).

README setup frame 1: 3.279 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    59.461 ms  35.676 Mcyc 100.0%
  pov_preserve_half      136.13 us   81.68 kcyc   0.2% x0.96875 141 us/call
  rs_draw_rings          23.984 ms  14.390 Mcyc  40.3%
    rs_ring_scan         23.272 ms  13.963 Mcyc  39.1% x54.625 426 us/call
  rs_timeline_step        60.88 us   36.53 kcyc   0.1% x1 61 us/call
  canvas_clear            89.17 us   53.50 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     35.181 ms  21.109 Mcyc  59.2% x1 35181 us/call
```

Wall min/avg/max = 3.279/59.460/87.570 ms. `rs_draw_rings` averages 23.984 ms/f; total render averages 24.279 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 737–768)

```
frame                    62.596 ms  37.557 Mcyc 100.0%
  pov_preserve_half      141.27 us   84.76 kcyc   0.2% x1 141 us/call
  rs_draw_rings          39.411 ms  23.646 Mcyc  63.0%
    rs_ring_scan         38.387 ms  23.032 Mcyc  61.3% x76 505 us/call
  rs_timeline_step        61.91 us   37.14 kcyc   0.1% x1 62 us/call
  canvas_clear            85.54 us   51.32 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     22.887 ms  13.732 Mcyc  36.6% x1 22887 us/call
```

Wall min/avg/max = 44.344/62.595/80.459 ms. `rs_draw_rings` averages 39.411 ms/f; total render averages 39.709 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.684 ms  37.610 Mcyc 100.0%
  pov_preserve_half      141.63 us   84.98 kcyc   0.2% x1 142 us/call
  rs_draw_rings          37.588 ms  22.553 Mcyc  60.0%
    rs_ring_scan         36.548 ms  21.929 Mcyc  58.3% x76 481 us/call
  rs_timeline_step        61.78 us   37.07 kcyc   0.1% x1 62 us/call
  canvas_clear            85.44 us   51.26 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     24.797 ms  14.878 Mcyc  39.6% x1 24797 us/call
```

Wall min/avg/max = 56.652/62.683/69.737 ms. `rs_draw_rings` averages 37.588 ms/f; total render averages 37.887 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 737–768; CPU shares use its separately snapshotted ISR interval of 2004771 us.

```
isr_wake       1154.8/f 0.375/1.645/20.586 us 3.03%
  isr_pack       144.3/f 5.613/7.883/10.945 us 1.81%
  isr_dma_submit 144.3/f 0.768/0.913/3.766 us 0.21%
```

- Mean pack 7.883 us/call versus submit 0.913 us/call (8.6×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.03%, corresponding to 60.606 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.388×; Peak live show render: 0.635×; Show end: 0.606×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `rs_draw_rings` — 63.0% of the frame, 39.411 ms/f.
2. `rs_ring_scan` — 61.3% of the frame, 38.387 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.141 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared Plot raster, column-walk and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
