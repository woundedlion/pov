# Comets on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\comets_ship.log`. Replaces `profile_comets_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the shared Plot raster, column-walk and filter-pipeline kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Comets 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 120 s / 1920 nominal frames at 16 fps; window 16. Raw capture 135 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Comets profile 123 16`; analyze only the 120 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:89760, data:149148, headers:8900; free:1783808 / RAM1: variables:313952, code:26952, padding:5816; free:177568 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 641–656 root cycles / 600 MHz = 992513.200 us; measured wall sum 992512 us, delta 1.21 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 09:29, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `3ae391de713860ba09dcbf5c6df83468b202b19ffa5b6459441da406443b8fa6`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 120 s (1920 nominal frames at 16 fps). Show-frame peak render is **29.869 ms**, spilled **0/1920**. Setup frame 1 is **0.869 ms** and remains separate. Live show peak is **29.869 ms**, spilled **0/1919**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1920. The raw log retains 2137 frame records; the 217 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 29.869 (12), spilled 🟢 0/1919 (0.00%).

README setup frame 1: 0.869 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    54.862 ms  32.917 Mcyc 100.0%
  pov_preserve_half      139.93 us   83.96 kcyc   0.3% x0.9375 149 us/call
  cm_draw_trail           1.763 ms   1.058 Mcyc   3.2% x1 1763 us/call
  cm_wipe_rebake           0.09 us    0.05 kcyc   0.0% x1 0 us/call
  cm_timeline_step        94.71 us   56.82 kcyc   0.2% x1 95 us/call
  canvas_clear            93.64 us   56.18 kcyc   0.2% x1 94 us/call
  canvas_buffer_wait     52.760 ms  31.656 Mcyc  96.2% x1 52760 us/call
```

Wall min/avg/max = 0.870/54.861/63.457 ms. `cm_draw_trail` averages 1.763 ms/f; total render averages 2.102 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 641–656)

```
frame                    62.032 ms  37.219 Mcyc 100.0%
  pov_preserve_half      138.09 us   82.85 kcyc   0.2% x1 138 us/call
  cm_draw_trail          19.751 ms  11.851 Mcyc  31.8% x1 19751 us/call
  cm_wipe_rebake          1.757 ms   1.054 Mcyc   2.8% x1 1757 us/call
  cm_timeline_step       141.05 us   84.63 kcyc   0.2% x1 141 us/call
  canvas_clear            84.77 us   50.86 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     40.149 ms  24.090 Mcyc  64.7% x1 40149 us/call
```

Wall min/avg/max = 45.610/62.032/77.352 ms. `cm_draw_trail` averages 19.751 ms/f; total render averages 21.883 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Costliest qualified preset (frames 513–528)

```
frame                    62.486 ms  37.491 Mcyc 100.0%
  pov_preserve_half      139.72 us   83.83 kcyc   0.2% x1 140 us/call
  cm_draw_trail          20.271 ms  12.163 Mcyc  32.4% x1 20271 us/call
  cm_wipe_rebake          1.754 ms   1.052 Mcyc   2.8% x1 1754 us/call
  cm_timeline_step       153.69 us   92.22 kcyc   0.2% x1 154 us/call
  canvas_clear            84.98 us   50.99 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     40.072 ms  24.043 Mcyc  64.1% x1 40072 us/call
```

Wall min/avg/max = 58.271/62.485/66.761 ms. `cm_draw_trail` averages 20.271 ms/f; total render averages 22.414 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 625–640)

```
frame                    62.000 ms  37.200 Mcyc 100.0%
  pov_preserve_half      139.90 us   83.94 kcyc   0.2% x1 140 us/call
  cm_draw_trail          19.714 ms  11.828 Mcyc  31.8% x1 19714 us/call
  cm_wipe_rebake           0.09 us    0.05 kcyc   0.0% x1 0 us/call
  cm_timeline_step       135.48 us   81.29 kcyc   0.2% x1 135 us/call
  canvas_clear            84.68 us   50.81 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     41.912 ms  25.147 Mcyc  67.6% x1 41912 us/call
```

Wall min/avg/max = 56.993/61.999/66.570 ms. `cm_draw_trail` averages 19.714 ms/f; total render averages 20.088 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1905–1920)

```
frame                    62.251 ms  37.350 Mcyc 100.0%
  pov_preserve_half      140.23 us   84.14 kcyc   0.2% x1 140 us/call
  cm_draw_trail          14.748 ms   8.849 Mcyc  23.7% x1 14748 us/call
  cm_wipe_rebake           0.09 us    0.05 kcyc   0.0% x1 0 us/call
  cm_timeline_step       137.32 us   82.39 kcyc   0.2% x1 137 us/call
  canvas_clear            84.65 us   50.79 kcyc   0.1% x1 85 us/call
  canvas_buffer_wait     47.125 ms  28.275 Mcyc  75.7% x1 47125 us/call
```

Wall min/avg/max = 59.837/62.250/64.304 ms. `cm_draw_trail` averages 14.748 ms/f; total render averages 15.126 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Held windows use fixed authored parameters while the scene keeps animating. Qualified windows exclude incoming parameter blends, fade-in and fade-out, and mixed-owner windows.

The show-scoped data covers 12 preset identities. Completing the scheduled show establishes coverage. Cadence buckets include all live frames, including in-show transitions. Qualified windows start after setup frame 1. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 4 | 9 | 513–528 | held | — | — | 20.271 | 22.414 | 16 | 26.006 | 0/160 |
| 10 | 9 | 1489–1504 | held | — | — | 19.931 | 20.295 | 16 | 26.414 | 0/160 |
| 12 | 9 | 1857–1872 | held | — | — | 19.847 | 20.208 | 16 | 22.891 | 0/160 |
| 5 | 9 | 641–656 | held | — | — | 19.751 | 21.883 | 16 | 29.869 | 0/160 |
| 2 | 9 | 289–304 | held | — | — | 19.439 | 19.796 | 16 | 24.835 | 0/160 |
| 8 | 9 | 1233–1248 | held | — | — | 19.220 | 19.589 | 16 | 22.088 | 0/160 |
| 9 | 9 | 1313–1328 | held | — | — | 18.746 | 20.898 | 16 | 25.253 | 0/160 |
| 11 | 9 | 1633–1648 | held | — | — | 17.776 | 19.920 | 16 | 26.775 | 0/160 |
| 6 | 9 | 801–816 | held | — | — | 14.774 | 16.852 | 16 | 20.416 | 0/160 |
| 3 | 9 | 433–448 | held | — | — | 14.757 | 15.099 | 16 | 18.595 | 0/160 |
| 7 | 9 | 1089–1104 | held | — | — | 12.412 | 12.745 | 16 | 15.579 | 0/160 |
| 1 | 8 | 113–128 | held | — | — | 12.128 | 12.451 | 16 | 14.636 | 0/159 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 641–656; CPU shares use its separately snapshotted ISR interval of 994181 us.

```
isr_wake       1145.4/f 0.431/1.458/13.240 us 2.69%
  isr_pack       143.2/f 5.506/6.545/10.228 us 1.50%
  isr_dma_submit 143.2/f 0.663/0.913/1.185 us 0.21%
```

- Mean pack 6.545 us/call versus submit 0.913 us/call (7.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 2.69%, corresponding to 60.819 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.034×; Peak live show render: 0.350×; Costliest qualified preset: 0.359×; Transition / non-hold regime: 0.321×; Show end: 0.242×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `cm_draw_trail` — 31.8% of the frame, 19.751 ms/f.
2. `cm_wipe_rebake` — 2.8% of the frame, 1.757 ms/f.
3. `cm_timeline_step` — 0.2% of the frame, 0.141 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the shared Plot raster, column-walk and filter-pipeline kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Comets`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
