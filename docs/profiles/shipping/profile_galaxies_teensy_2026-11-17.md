# Galaxies on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\galaxies_ship.log`. Replaces `profile_galaxies_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared scan-shapes, spherical-math and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Galaxies 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 60 s / 960 nominal frames at 16 fps; window 16. Raw capture 63 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Galaxies profile 63 16`; analyze only the 60 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:92080, data:149436, headers:8340; free:1781760 / RAM1: variables:313952, code:23672, padding:9096; free:177568 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 865–880 root cycles / 600 MHz = 997549.000 us; measured wall sum 997548 us, delta 1.00 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:28, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `f8fb23302c9c20249b95d6dc957a5a4eadf82dc9bd61262870db04b8beef50c6`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 60 s (960 nominal frames at 16 fps). Show-frame peak render is **43.803 ms**, spilled **0/960**. Setup frame 1 is **1.477 ms** and remains separate. Live show peak is **43.803 ms**, spilled **0/959**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–960. The raw log retains 984 frame records; the 24 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 43.803, spilled 🟢 0/959 (0.00%).

README setup frame 1: 1.477 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    54.803 ms  32.882 Mcyc 100.0%
  pov_preserve_half      139.43 us   83.66 kcyc   0.3% x0.9375 149 us/call
  gx_draw_particles      773.86 us  464.32 kcyc   1.4% x1 774 us/call
  gx_particle_step       400.38 us  240.23 kcyc   0.7% x1 400 us/call
  gx_timeline_step        49.94 us   29.97 kcyc   0.1% x1 50 us/call
  canvas_clear            93.81 us   56.29 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     53.327 ms  31.996 Mcyc  97.3% x1 53327 us/call
```

Wall min/avg/max = 1.219/54.803/62.681 ms. `gx_particle_step` averages 0.400 ms/f; total render averages 1.477 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 865–880)

```
frame                    62.347 ms  37.408 Mcyc 100.0%
  pov_preserve_half      143.94 us   86.36 kcyc   0.2% x1 144 us/call
  gx_draw_particles      15.935 ms   9.561 Mcyc  25.6% x1 15935 us/call
  gx_particle_step       26.884 ms  16.130 Mcyc  43.1% x1 26884 us/call
  gx_timeline_step        49.38 us   29.63 kcyc   0.1% x1 49 us/call
  canvas_clear            85.58 us   51.35 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     19.228 ms  11.537 Mcyc  30.8% x1 19228 us/call
```

Wall min/avg/max = 61.136/62.346/63.702 ms. `gx_particle_step` averages 26.884 ms/f; total render averages 43.119 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 945–960)

```
frame                    62.388 ms  37.433 Mcyc 100.0%
  pov_preserve_half      144.45 us   86.67 kcyc   0.2% x1 144 us/call
  gx_draw_particles      15.144 ms   9.087 Mcyc  24.3% x1 15144 us/call
  gx_particle_step       26.834 ms  16.100 Mcyc  43.0% x1 26834 us/call
  gx_timeline_step        47.87 us   28.72 kcyc   0.1% x1 48 us/call
  canvas_clear            85.71 us   51.43 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     20.113 ms  12.068 Mcyc  32.2% x1 20113 us/call
```

Wall min/avg/max = 60.496/62.388/64.311 ms. `gx_particle_step` averages 26.834 ms/f; total render averages 42.275 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 865–880; CPU shares use its separately snapshotted ISR interval of 999506 us.

```
isr_wake       1151.5/f 0.296/1.561/17.606 us 2.88%
  isr_pack       143.9/f 5.730/7.238/10.155 us 1.66%
  isr_dma_submit 143.9/f 0.623/0.910/0.975 us 0.20%
```

- Mean pack 7.238 us/call versus submit 0.910 us/call (8.0×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.88%, corresponding to 60.700 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.024×; Peak live show render: 0.690×; Show end: 0.676×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `gx_particle_step` — 43.1% of the frame, 26.884 ms/f.
2. `gx_draw_particles` — 25.6% of the frame, 15.935 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.144 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared scan-shapes, spherical-math and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Galaxies`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
