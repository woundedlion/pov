# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\fishbowl_ship.log`. Replaces `profile_fishbowl_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared scan-shapes, spherical-math and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 32. Raw capture 70 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Fishbowl profile 64 32`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:139912, data:149284, headers:8788; free:1733632 / RAM1: variables:314048, code:37864, padding:27672; free:144704 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 193–224 root cycles / 600 MHz = 1996230.997 us; measured wall sum 1996231 us, delta 0.00 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:34, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `1761b23d9cd41a25bb8aaad61fc0a2df2fef469f302a0577168e7518e4a3fd55`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **28.256 ms**, spilled **0/960**. Setup frame 1 is **1.681 ms** and remains separate. Live show peak is **28.256 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 1096 frame records; the 136 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 28.256, spilled 🟢 0/959 (0.00%).

README setup frame 1: 1.681 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    58.782 ms  35.269 Mcyc 100.0%
  pov_preserve_half      140.89 us   84.53 kcyc   0.2% x0.96875 145 us/call
  fish_multiline_draw     1.852 ms   1.111 Mcyc   3.2% x1 1852 us/call
  fish_build_vertices     1.844 ms   1.106 Mcyc   3.1% x1 1844 us/call
  fish_noise_prepare       4.16 us    2.50 kcyc   0.0% x1 4 us/call
  fish_timeline_step     186.18 us  111.71 kcyc   0.3%
    animation_noise_clock    0.18 us    0.11 kcyc   0.0% x1 0 us/call
  canvas_clear            94.42 us   56.65 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     54.492 ms  32.695 Mcyc  92.7% x1 54492 us/call
```

Wall min/avg/max = 1.354/58.782/64.013 ms. `fish_build_vertices` averages 1.844 ms/f; total render averages 4.290 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 193–224)

```
frame                    62.382 ms  37.429 Mcyc 100.0%
  pov_preserve_half      137.66 us   82.60 kcyc   0.2% x1 138 us/call
  fish_multiline_draw    10.533 ms   6.320 Mcyc  16.9% x1 10533 us/call
  fish_build_vertices    12.538 ms   7.523 Mcyc  20.1% x1 12538 us/call
  fish_noise_prepare       4.57 us    2.74 kcyc   0.0% x1 5 us/call
  fish_timeline_step     202.27 us  121.36 kcyc   0.3%
    animation_noise_clock    0.18 us    0.11 kcyc   0.0% x1 0 us/call
  canvas_clear            91.94 us   55.17 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     38.234 ms  22.941 Mcyc  61.3% x1 38234 us/call
```

Wall min/avg/max = 57.537/62.382/67.823 ms. `fish_build_vertices` averages 12.538 ms/f; total render averages 24.148 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 929–960)

```
frame                    62.358 ms  37.415 Mcyc 100.0%
  pov_preserve_half      140.81 us   84.48 kcyc   0.2% x1 141 us/call
  fish_multiline_draw     9.260 ms   5.556 Mcyc  14.8% x1 9260 us/call
  fish_build_vertices    12.563 ms   7.538 Mcyc  20.1% x1 12563 us/call
  fish_noise_prepare       4.90 us    2.94 kcyc   0.0% x1 5 us/call
  fish_timeline_step     208.80 us  125.28 kcyc   0.3%
    animation_noise_clock    0.18 us    0.11 kcyc   0.0% x1 0 us/call
  canvas_clear            92.12 us   55.27 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait     39.440 ms  23.664 Mcyc  63.2% x1 39440 us/call
```

Wall min/avg/max = 59.023/62.357/66.446 ms. `fish_build_vertices` averages 12.563 ms/f; total render averages 22.917 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 193–224; CPU shares use its separately snapshotted ISR interval of 2001087 us.

```
isr_wake       1152.7/f 0.270/1.376/20.720 us 2.53%
  isr_pack       144.1/f 5.540/6.240/8.951 us 1.43%
  isr_dma_submit 144.1/f 0.606/0.913/8.710 us 0.21%
```

- Mean pack 6.240 us/call versus submit 0.913 us/call (6.8×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.53%, corresponding to 60.919 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.069×; Peak live show render: 0.386×; Show end: 0.367×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `fish_build_vertices` — 20.1% of the frame, 12.538 ms/f.
2. `fish_multiline_draw` — 16.9% of the frame, 10.533 ms/f.
3. `fish_timeline_step` — 0.3% of the frame, 0.202 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared scan-shapes, spherical-math and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Fishbowl`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
