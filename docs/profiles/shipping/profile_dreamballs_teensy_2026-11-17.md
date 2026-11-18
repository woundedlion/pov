# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\dreamballs_ship.log`. Replaces `profile_dreamballs_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared Plot raster, column-walk and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | DreamBalls 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 200 s / 3200 nominal frames at 16 fps; window 16. Raw capture 215 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh DreamBalls profile 203 16`; analyze only the 200 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:154360, data:190212, headers:8708; free:1678336 / RAM1: variables:314144, code:36408, padding:29128; free:144608 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 2737–2752 root cycles / 600 MHz = 1004752.098 us; measured wall sum 1004752 us, delta 0.10 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:19, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `f5d364140086744d2acab22d91ea1892bf6ad6741882bf0dcc040dfe4fd67ff3`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 200 s (3200 nominal frames at 16 fps). Show-frame peak render is **51.363 ms**, spilled **0/3200**. Setup frame 1 is **49.416 ms** and remains separate. Live show peak is **51.363 ms**, spilled **0/3199**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–3200. The raw log retains 3416 frame records; the 216 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 51.363 (10), spilled 🟢 0/3199 (0.00%).

README setup frame 1: 49.416 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    60.205 ms  36.123 Mcyc 100.0%
  pov_preserve_half      136.45 us   81.87 kcyc   0.2% x0.9375 146 us/call
  db_timeline_step       38.416 ms  23.050 Mcyc  63.8%
    db_draw              38.345 ms  23.007 Mcyc  63.7%
      db_draw_scene      38.329 ms  22.997 Mcyc  63.7%
        db_mesh_plot     37.786 ms  22.672 Mcyc  62.8% x18 2099 us/call
        db_displace      463.73 us  278.24 kcyc   0.8% x18 26 us/call
  canvas_clear            94.04 us   56.42 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     21.546 ms  12.927 Mcyc  35.8% x1 21546 us/call
```

Wall min/avg/max = 36.331/60.204/66.733 ms. `db_timeline_step` averages 38.416 ms/f; total render averages 38.659 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 2737–2752)

```
frame                    62.797 ms  37.678 Mcyc 100.0%
  pov_preserve_half      143.16 us   85.90 kcyc   0.2% x1 143 us/call
  db_timeline_step       47.915 ms  28.749 Mcyc  76.3%
    db_draw              47.783 ms  28.670 Mcyc  76.1%
      db_draw_scene      47.767 ms  28.660 Mcyc  76.1%
        db_mesh_plot     46.383 ms  27.830 Mcyc  73.9% x6 7730 us/call
        db_displace       1.129 ms  677.56 kcyc   1.8% x6 188 us/call
  canvas_clear            85.90 us   51.54 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     14.636 ms   8.782 Mcyc  23.3% x1 14636 us/call
```

Wall min/avg/max = 59.367/62.797/65.375 ms. `db_timeline_step` averages 47.915 ms/f; total render averages 48.161 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 2769–2784)

```
frame                    62.340 ms  37.404 Mcyc 100.0%
  pov_preserve_half      142.25 us   85.35 kcyc   0.2% x1 142 us/call
  db_timeline_step       48.244 ms  28.946 Mcyc  77.4%
    db_draw              48.105 ms  28.863 Mcyc  77.2%
      db_draw_scene      48.089 ms  28.853 Mcyc  77.1%
        db_mesh_plot     46.691 ms  28.014 Mcyc  74.9% x6 7782 us/call
        db_displace       1.130 ms  678.21 kcyc   1.8% x6 188 us/call
  canvas_clear            86.02 us   51.61 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     13.847 ms   8.308 Mcyc  22.2% x1 13847 us/call
```

Wall min/avg/max = 58.775/62.340/65.373 ms. `db_timeline_step` averages 48.244 ms/f; total render averages 48.493 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 2577–2592)

```
frame                    62.487 ms  37.492 Mcyc 100.0%
  pov_preserve_half      143.68 us   86.21 kcyc   0.2% x1 144 us/call
  db_timeline_step       45.767 ms  27.460 Mcyc  73.2%
    db_draw              45.652 ms  27.391 Mcyc  73.1%
      db_draw_scene      45.637 ms  27.382 Mcyc  73.0%
        db_mesh_plot     44.268 ms  26.561 Mcyc  70.8% x6 7378 us/call
        db_displace       1.123 ms  673.60 kcyc   1.8% x6 187 us/call
  canvas_clear            85.42 us   51.25 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     16.477 ms   9.886 Mcyc  26.4% x1 16477 us/call
```

Wall min/avg/max = 57.400/62.486/67.431 ms. `db_timeline_step` averages 45.767 ms/f; total render averages 46.009 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 3185–3200)

```
frame                    60.422 ms  36.253 Mcyc 100.0%
  pov_preserve_half      144.05 us   86.43 kcyc   0.2% x1 144 us/call
  db_timeline_step       31.050 ms  18.630 Mcyc  51.4%
    db_draw              30.933 ms  18.560 Mcyc  51.2%
      db_draw_scene      30.918 ms  18.551 Mcyc  51.2%
        db_mesh_plot     30.252 ms  18.151 Mcyc  50.1% x5.625 5378 us/call
        db_displace      531.22 us  318.73 kcyc   0.9% x5.625 94 us/call
  canvas_clear            84.79 us   50.87 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     29.131 ms  17.478 Mcyc  48.2% x1 29131 us/call
```

Wall min/avg/max = 29.522/60.421/63.901 ms. `db_timeline_step` averages 31.050 ms/f; total render averages 31.291 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held windows use fixed authored parameters while the scene keeps animating. Qualified windows exclude incoming parameter blends, fade-in and fade-out, and mixed-owner windows. Presets change by snap. Sprites live for 320 frames with no overlap; the held windows exclude their independent 32-frame fade-in and fade-out, including the initial fade-in.

The show-scoped data covers 10 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 9 | 15 | 2769–2784 | held | — | — | 48.244 | 48.493 | 16 | 51.363 | 0/320 |
| 1 | 15 | 209–224 | held | — | — | 44.326 | 44.574 | 16 | 47.619 | 0/319 |
| 10 | 15 | 3089–3104 | held | — | — | 35.986 | 36.234 | 16 | 39.932 | 0/320 |
| 8 | 15 | 2417–2432 | held | — | — | 35.487 | 35.737 | 16 | 38.354 | 0/320 |
| 4 | 15 | 1137–1152 | held | — | — | 27.465 | 27.712 | 16 | 31.519 | 0/320 |
| 3 | 15 | 849–864 | held | — | — | 25.782 | 26.031 | 16 | 28.556 | 0/320 |
| 2 | 15 | 353–368 | held | — | — | 24.809 | 25.053 | 16 | 28.084 | 0/320 |
| 7 | 15 | 1953–1968 | held | — | — | 23.432 | 23.677 | 16 | 25.371 | 0/320 |
| 6 | 15 | 1633–1648 | held | — | — | 17.213 | 17.463 | 16 | 19.965 | 0/320 |
| 5 | 15 | 1473–1488 | held | — | — | 16.583 | 16.827 | 16 | 18.158 | 0/320 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 2737–2752; CPU shares use its separately snapshotted ISR interval of 1006438 us.

```
isr_wake       1159.4/f 0.390/1.606/17.380 us 2.96%
  isr_pack       144.9/f 6.180/7.393/12.148 us 1.70%
  isr_dma_submit 144.9/f 0.665/0.923/10.003 us 0.21%
```

- Mean pack 7.393 us/call versus submit 0.923 us/call (8.0×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.96%, corresponding to 60.650 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.619×; Peak live show render: 0.771×; Costliest qualified preset: 0.776×; Transition / non-hold regime: 0.736×; Show end: 0.501×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `db_timeline_step` — 76.3% of the frame, 47.915 ms/f.
2. `db_draw` — 76.1% of the frame, 47.783 ms/f.
3. `db_draw_scene` — 76.1% of the frame, 47.767 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared Plot raster, column-walk and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
