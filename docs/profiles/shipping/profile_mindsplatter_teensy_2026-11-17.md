# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\mindsplatter_ship.log`. Replaces `profile_mindsplatter_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared Plot raster, column-walk and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 80 s / 1280 nominal frames at 16 fps; window 16. Raw capture 98 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh MindSplatter profile 83 16`; analyze only the 80 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:156944, data:201916, headers:8732; free:1664024 / RAM1: variables:314176, code:45384, padding:20152; free:144576 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 993–1008 root cycles / 600 MHz = 1001762.160 us; measured wall sum 1001762 us, delta 0.16 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:03, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `533aa2bdd8d6827a168ceeec8349e82b21d979f59e8ac117d462ffd90749d9a2`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 80 s (1280 nominal frames at 16 fps). Show-frame peak render is **50.810 ms**, spilled **0/1280**. Setup frame 1 is **0.372 ms** and remains separate. Live show peak is **50.810 ms**, spilled **0/1279**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1280. The raw log retains 1522 frame records; the 242 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 50.810 (8), spilled 🟢 0/1279 (0.00%).

README setup frame 1: 0.372 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    54.764 ms  32.858 Mcyc 100.0%
  pov_preserve_half      137.19 us   82.31 kcyc   0.3% x0.9375 146 us/call
  msp_draw_particles     528.12 us  316.87 kcyc   1.0%
    msp_particle_scan      1.19 us    0.71 kcyc   0.0% x0.0625 19 us/call
  msp_particle_step      124.91 us   74.94 kcyc   0.2% x1 125 us/call
  msp_timeline_step       57.61 us   34.56 kcyc   0.1% x1 58 us/call
  canvas_clear            92.05 us   55.23 kcyc   0.2% x1 92 us/call
  canvas_buffer_wait     53.810 ms  32.286 Mcyc  98.3% x1 53810 us/call
```

Wall min/avg/max = 0.374/54.763/62.598 ms. `msp_draw_particles` averages 0.528 ms/f; total render averages 0.954 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 993–1008)

```
frame                    62.610 ms  37.566 Mcyc 100.0%
  pov_preserve_half      136.99 us   82.19 kcyc   0.2% x1 137 us/call
  msp_draw_particles     42.669 ms  25.601 Mcyc  68.1% x1 42669 us/call
  msp_particle_step       2.956 ms   1.774 Mcyc   4.7% x1 2956 us/call
  msp_timeline_step      124.95 us   74.97 kcyc   0.2% x1 125 us/call
  canvas_clear            89.48 us   53.69 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     16.594 ms   9.957 Mcyc  26.5% x1 16594 us/call
```

Wall min/avg/max = 59.172/62.610/67.345 ms. `msp_draw_particles` averages 42.669 ms/f; total render averages 46.016 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 1009–1024)

```
frame                    61.694 ms  37.016 Mcyc 100.0%
  pov_preserve_half      137.62 us   82.57 kcyc   0.2% x1 138 us/call
  msp_draw_particles     35.828 ms  21.497 Mcyc  58.1% x1 35828 us/call
  msp_particle_step       2.668 ms   1.601 Mcyc   4.3% x1 2668 us/call
  msp_timeline_step      108.42 us   65.05 kcyc   0.2% x1 108 us/call
  canvas_clear            88.18 us   52.91 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     22.823 ms  13.694 Mcyc  37.0% x1 22823 us/call
```

Wall min/avg/max = 53.377/61.693/69.173 ms. `msp_draw_particles` averages 35.828 ms/f; total render averages 38.871 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1265–1280)

```
frame                    62.560 ms  37.536 Mcyc 100.0%
  pov_preserve_half      140.62 us   84.37 kcyc   0.2% x1 141 us/call
  msp_draw_particles     26.826 ms  16.095 Mcyc  42.9% x1 26826 us/call
  msp_particle_step       2.075 ms   1.245 Mcyc   3.3% x1 2075 us/call
  msp_timeline_step      114.05 us   68.43 kcyc   0.2% x1 114 us/call
  canvas_clear            89.08 us   53.45 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     33.267 ms  19.960 Mcyc  53.2% x1 33267 us/call
```

Wall min/avg/max = 60.281/62.559/65.460 ms. `msp_draw_particles` averages 26.826 ms/f; total render averages 29.293 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held windows use fixed authored parameters while the scene keeps animating. Qualified windows exclude incoming parameter blends, fade-in and fade-out, and mixed-owner windows.

The show-scoped data covers 8 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 7 | 6 | 1009–1024 | held | — | — | 35.828 | 38.871 | 16 | 50.810 | 0/159 |
| 6 | 6 | 929–944 | held | — | — | 35.350 | 43.453 | 16 | 48.917 | 0/159 |
| 4 | 6 | 529–544 | held | — | — | 29.590 | 32.980 | 16 | 49.826 | 0/159 |
| 2 | 6 | 209–224 | held | — | — | 28.724 | 31.274 | 16 | 49.078 | 0/159 |
| 8 | 6 | 1233–1248 | held | — | — | 25.283 | 27.750 | 16 | 31.188 | 0/159 |
| 3 | 6 | 369–384 | held | — | — | 23.971 | 26.375 | 16 | 30.138 | 0/159 |
| 1 | 8 | 129–144 | held | — | — | 17.824 | 20.118 | 16 | 30.609 | 0/166 |
| 5 | 6 | 689–704 | held | — | — | 15.465 | 17.133 | 16 | 30.223 | 0/159 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 993–1008; CPU shares use its separately snapshotted ISR interval of 1003203 us.

```
isr_wake       1155.8/f 0.276/1.620/21.153 us 2.98%
  isr_pack       144.4/f 5.740/7.633/10.533 us 1.75%
  isr_dma_submit 144.4/f 0.655/0.915/6.365 us 0.21%
```

- Mean pack 7.633 us/call versus submit 0.915 us/call (8.3×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.98%, corresponding to 60.637 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.015×; Peak live show render: 0.736×; Costliest qualified preset: 0.622×; Show end: 0.469×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `msp_draw_particles` — 68.1% of the frame, 42.669 ms/f.
2. `msp_particle_step` — 4.7% of the frame, 2.956 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.137 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared Plot raster, column-walk and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
