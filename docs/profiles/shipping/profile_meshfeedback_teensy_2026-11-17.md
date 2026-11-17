# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\meshfeedback_ship.log`. Replaces `profile_meshfeedback_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 181 s / 2896 nominal frames at 16 fps; window 16. Raw capture 196 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh MeshFeedback profile 184 16`; analyze only the 181 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:195344, data:190696, headers:8200; free:1637376 / RAM1: variables:314176, code:45976, padding:19560; free:144576 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 2641–2656 root cycles / 600 MHz = 1001521.397 us; measured wall sum 1001521 us, delta 0.40 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:23, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `31eded98152be49273c6b1b973b8ba7f048a034b477fc915abf759f263e564f8`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 181 s (2896 nominal frames at 16 fps). Show-frame peak render is **55.672 ms**, spilled **0/2896**. Setup frame 1 is **6.470 ms** and remains separate. Live show peak is **55.672 ms**, spilled **0/2895**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the full canvas, containing 41,472 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–2896. The raw log retains 3112 frame records; the 216 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 55.672 (12), spilled 🟢 0/2895 (0.00%).

README setup frame 1: 6.470 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    57.893 ms  34.736 Mcyc 100.0%
  feedback_forward        1.566 ms  939.79 kcyc   2.7% x0.9375 1671 us/call
  mf_mesh_draw            6.774 ms   4.064 Mcyc  11.7% x1 6774 us/call
  mf_feedback_flush      32.735 ms  19.641 Mcyc  56.5%
    feedback_composite   24.119 ms  14.472 Mcyc  41.7% x0.9375 25727 us/call
    feedback_populate     7.197 ms   4.318 Mcyc  12.4% x0.9375 7676 us/call
    feedback_litscan      64.30 us   38.58 kcyc   0.1% x1 64 us/call
  mf_timeline_step        95.33 us   57.20 kcyc   0.2%
    animation_noise_clock    0.18 us    0.11 kcyc   0.0% x1 0 us/call
  mf_apply_params          2.58 us    1.55 kcyc   0.0% x1 3 us/call
  canvas_buffer_wait     16.652 ms   9.991 Mcyc  28.8% x1 16652 us/call
```

Wall min/avg/max = 6.473/57.893/65.716 ms. `mf_feedback_flush` averages 32.735 ms/f; total render averages 41.241 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 2641–2656)

```
frame                    62.595 ms  37.557 Mcyc 100.0%
  feedback_forward       774.46 us  464.68 kcyc   1.2% x1 774 us/call
  mf_mesh_draw            9.435 ms   5.661 Mcyc  15.1% x1 9435 us/call
  mf_feedback_flush      34.116 ms  20.469 Mcyc  54.5%
    feedback_composite   33.136 ms  19.881 Mcyc  52.9% x1 33136 us/call
    feedback_populate    479.96 us  287.98 kcyc   0.8% x1 480 us/call
    feedback_litscan      13.23 us    7.94 kcyc   0.0% x1 13 us/call
  mf_timeline_step        93.69 us   56.21 kcyc   0.1% x1 94 us/call
  mf_apply_params         51.35 us   30.81 kcyc   0.1% x1 51 us/call
  canvas_buffer_wait     18.052 ms  10.831 Mcyc  28.8% x1 18052 us/call
```

Wall min/avg/max = 53.930/62.595/75.065 ms. `mf_feedback_flush` averages 34.116 ms/f; total render averages 44.543 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 513–528)

```
frame                    62.404 ms  37.442 Mcyc 100.0%
  feedback_forward        1.110 ms  666.10 kcyc   1.8% x1 1110 us/call
  mf_mesh_draw            7.730 ms   4.638 Mcyc  12.4% x1 7730 us/call
  mf_feedback_flush      42.598 ms  25.559 Mcyc  68.3%
    feedback_composite   33.389 ms  20.033 Mcyc  53.5% x1 33389 us/call
    feedback_populate     7.750 ms   4.650 Mcyc  12.4% x1 7750 us/call
    feedback_litscan      12.29 us    7.37 kcyc   0.0% x1 12 us/call
  mf_timeline_step        96.28 us   57.77 kcyc   0.2%
    animation_noise_clock    0.18 us    0.11 kcyc   0.0% x1 0 us/call
  mf_apply_params          1.77 us    1.06 kcyc   0.0% x1 2 us/call
  canvas_buffer_wait     10.802 ms   6.481 Mcyc  17.3% x1 10802 us/call
```

Wall min/avg/max = 61.882/62.403/62.976 ms. `mf_feedback_flush` averages 42.598 ms/f; total render averages 51.602 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 1201–1216)

```
frame                    62.381 ms  37.429 Mcyc 100.0%
  feedback_forward        1.356 ms  813.36 kcyc   2.2% x1 1356 us/call
  mf_mesh_draw            7.925 ms   4.755 Mcyc  12.7% x1 7925 us/call
  mf_feedback_flush      41.682 ms  25.009 Mcyc  66.8%
    feedback_composite   32.523 ms  19.514 Mcyc  52.1% x1 32523 us/call
    feedback_populate     7.708 ms   4.625 Mcyc  12.4% x1 7708 us/call
    feedback_litscan       9.97 us    5.98 kcyc   0.0% x1 10 us/call
  mf_timeline_step        97.56 us   58.54 kcyc   0.2%
    animation_noise_clock    0.18 us    0.11 kcyc   0.0% x1 0 us/call
  mf_apply_params         19.43 us   11.66 kcyc   0.0% x1 19 us/call
  canvas_buffer_wait     11.233 ms   6.740 Mcyc  18.0% x1 11233 us/call
```

Wall min/avg/max = 60.576/62.381/65.626 ms. `mf_feedback_flush` averages 41.682 ms/f; total render averages 51.149 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 2881–2896)

```
frame                    62.543 ms  37.526 Mcyc 100.0%
  feedback_forward        1.424 ms  854.25 kcyc   2.3% x1 1424 us/call
  mf_mesh_draw           10.233 ms   6.140 Mcyc  16.4% x1 10233 us/call
  mf_feedback_flush      36.447 ms  21.868 Mcyc  58.3%
    feedback_composite   33.297 ms  19.978 Mcyc  53.2% x1 33297 us/call
    feedback_populate     2.422 ms   1.453 Mcyc   3.9% x1 2422 us/call
    feedback_litscan       6.74 us    4.04 kcyc   0.0% x1 7 us/call
  mf_timeline_step        96.68 us   58.01 kcyc   0.2%
    animation_noise_clock    0.05 us    0.03 kcyc   0.0% x0.3125 0 us/call
  mf_apply_params          9.19 us    5.51 kcyc   0.0% x1 9 us/call
  canvas_buffer_wait     14.262 ms   8.557 Mcyc  22.8% x1 14262 us/call
```

Wall min/avg/max = 60.356/62.543/68.169 ms. `mf_feedback_flush` averages 36.447 ms/f; total render averages 48.281 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held means fixed authored parameters, not motionless pixels. Exclude complete incoming parameter blends and both halves of fades; cadence buckets retain all transitions. Never substitute mixed-owner or transition windows for absent clean holds.

The show-scoped data covers 12 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 3 · EnergeticFire | 14 | 513–528 | held | — | — | 42.598 | 51.602 | 16 | 52.703 | 0/241 |
| 9 · Miasma | 14 | 2129–2144 | held | — | — | 42.492 | 52.630 | 16 | 53.569 | 0/241 |
| 5 · SlowDust | 14 | 1137–1152 | held | — | — | 42.298 | 51.109 | 16 | 51.849 | 0/241 |
| 4 · Smoke | 14 | 737–752 | held | — | — | 41.476 | 45.059 | 16 | 45.809 | 0/241 |
| 2 · SlowFire | 15 | 257–272 | held | — | — | 41.380 | 47.254 | 16 | 48.159 | 0/241 |
| 6 · WavyTrails | 14 | 1297–1312 | held | — | — | 41.002 | 50.735 | 16 | 53.319 | 0/241 |
| 8 · MeltingLo | 14 | 1809–1824 | held | — | — | 40.761 | 50.924 | 16 | 53.015 | 0/241 |
| 7 · MeltingHi | 14 | 1569–1584 | held | — | — | 38.561 | 44.174 | 16 | 48.525 | 0/241 |
| 1 · ArcingLightning | 15 | 17–32 | held | — | — | 37.214 | 46.547 | 16 | 52.668 | 0/244 |
| 10 · LooseWormhole | 14 | 2369–2384 | held | — | — | 34.900 | 42.313 | 16 | 50.662 | 0/241 |
| 12 · WigglingWormhole | 14 | 2865–2880 | held | — | — | 34.606 | 47.563 | 16 | 55.672 | 0/241 |
| 11 · TightWormhole | 14 | 2513–2528 | held | — | — | 34.482 | 44.665 | 16 | 52.621 | 0/241 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The full canvas has 41,472 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 2641–2656; CPU shares use its separately snapshotted ISR interval of 1003722 us.

```
isr_wake       1156.4/f 0.285/1.781/20.101 us 3.28%
  isr_pack       144.5/f 7.525/8.911/11.280 us 2.05%
  isr_dma_submit 144.5/f 0.758/0.916/3.786 us 0.21%
```

- Mean pack 8.911 us/call versus submit 0.916 us/call (9.7×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.28%, corresponding to 60.450 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.660×; Peak live show render: 0.713×; Costliest qualified preset: 0.826×; Transition / non-hold regime: 0.818×; Show end: 0.772×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `mf_feedback_flush` — 54.5% of the frame, 34.116 ms/f.
2. `feedback_composite` — 52.9% of the frame, 33.136 ms/f.
3. `mf_mesh_draw` — 15.1% of the frame, 9.435 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: `core/render/scan/mesh.h`, `face_raster.h` and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
