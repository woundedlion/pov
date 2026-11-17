# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\gsreactiondiffusion_ship.log`. Replaces `profile_gsreactiondiffusion_teensy_2026-11-05.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the effect-local physics/shading functions and shared scan/filter kernels |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 30 s / 480 nominal frames at 16 fps; window 32. Raw capture 130 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh GSReactionDiffusion profile 34 32`; analyze only the 30 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:84336, data:249032, headers:8648; free:1689600 / RAM1: variables:314240, code:17688, padding:15080; free:177280 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 385–416 root cycles / 600 MHz = 1996433.438 us; measured wall sum 1996436 us, delta 1.28 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:43, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `85436776265234631a6ec6373df7dbb0bad9de79713bdb3816f4251de0d1a0fd`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 30 s (480 nominal frames at 16 fps). Show-frame peak render is **47.331 ms**, spilled **0/480**. Setup frame 1 is **24.651 ms** and remains separate. Live show peak is **47.331 ms**, spilled **0/479**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–480. The raw log retains 2056 frame records; the 1576 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 47.331, spilled 🟢 0/479 (0.00%).

README setup frame 1: 24.651 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–32)

```
frame                    60.124 ms  36.075 Mcyc 100.0%
  pov_preserve_half      133.78 us   80.27 kcyc   0.2% x0.96875 138 us/call
  grd_render             24.269 ms  14.561 Mcyc  40.4%
    grd_shader_draw       9.540 ms   5.724 Mcyc  15.9% x6 1590 us/call
    grd_cull_flags       717.59 us  430.56 kcyc   1.2% x6 120 us/call
    grd_publish           1.950 ms   1.170 Mcyc   3.2% x6 325 us/call
    grd_physics           5.147 ms   3.088 Mcyc   8.6% x6 858 us/call
    grd_pigment           3.063 ms   1.838 Mcyc   5.1% x1 3063 us/call
    rd_lattice_basis     105.27 us   63.16 kcyc   0.2% x1 105 us/call
    grd_color_noise       2.980 ms   1.788 Mcyc   5.0% x1 2980 us/call
  rd_timeline_step        48.60 us   29.16 kcyc   0.1% x1 49 us/call
  canvas_clear            93.21 us   55.93 kcyc   0.2% x1 93 us/call
  canvas_buffer_wait     35.566 ms  21.340 Mcyc  59.2% x1 35566 us/call
grd_color_palette 430.63 us 258.38 kcyc 0.7% x0.03125 13780 us/call
grd_seed_reaction 442.10 us 265.26 kcyc 0.7% x0.03125 14147 us/call
```

Wall min/avg/max = 21.024/60.124/64.807 ms. `grd_render` averages 24.269 ms/f; total render averages 24.558 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 385–416)

```
frame                    62.389 ms  37.433 Mcyc 100.0%
  pov_preserve_half      135.14 us   81.09 kcyc   0.2% x1 135 us/call
  grd_render             46.321 ms  27.793 Mcyc  74.2%
    grd_shader_draw      26.978 ms  16.187 Mcyc  43.2% x6 4496 us/call
    grd_cull_flags       723.51 us  434.11 kcyc   1.2% x6 121 us/call
    grd_publish           4.646 ms   2.788 Mcyc   7.4% x6 774 us/call
    grd_physics           5.169 ms   3.101 Mcyc   8.3% x6 862 us/call
    grd_pigment           4.218 ms   2.531 Mcyc   6.8% x1 4218 us/call
    rd_lattice_basis     103.76 us   62.26 kcyc   0.2% x1 104 us/call
    grd_color_noise       3.000 ms   1.800 Mcyc   4.8% x1 3000 us/call
  rd_timeline_step        46.21 us   27.72 kcyc   0.1% x1 46 us/call
  canvas_clear            91.15 us   54.69 kcyc   0.1% x1 91 us/call
  canvas_buffer_wait     15.779 ms   9.468 Mcyc  25.3% x1 15779 us/call
```

Wall min/avg/max = 61.545/62.388/63.252 ms. `grd_render` averages 46.321 ms/f; total render averages 46.609 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 449–480)

```
frame                    62.048 ms  37.229 Mcyc 100.0%
  pov_preserve_half      135.90 us   81.54 kcyc   0.2% x1 136 us/call
  grd_render             27.717 ms  16.630 Mcyc  44.7%
    grd_shader_draw      10.912 ms   6.547 Mcyc  17.6% x5.65625 1929 us/call
    grd_cull_flags       678.01 us  406.80 kcyc   1.1% x5.65625 120 us/call
    grd_publish           2.966 ms   1.780 Mcyc   4.8% x5.625 527 us/call
    grd_physics           4.989 ms   2.994 Mcyc   8.0% x5.8125 858 us/call
    grd_pigment           3.025 ms   1.815 Mcyc   4.9% x0.96875 3123 us/call
    rd_lattice_basis     104.46 us   62.68 kcyc   0.2% x1 104 us/call
    grd_color_noise       2.989 ms   1.793 Mcyc   4.8% x1 2989 us/call
  rd_timeline_step        56.19 us   33.71 kcyc   0.1% x1 56 us/call
  canvas_clear            87.90 us   52.74 kcyc   0.1% x1 88 us/call
  canvas_buffer_wait     34.036 ms  20.422 Mcyc  54.9% x1 34036 us/call
grd_color_palette        472.85 us  283.71 kcyc   0.8% x0.1875 2522 us/call
grd_seed_reaction        453.91 us  272.35 kcyc   0.7% x0.0625 7263 us/call
```

Wall min/avg/max = 58.816/62.048/64.777 ms. `grd_render` averages 27.717 ms/f; total render averages 28.012 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 385–416; CPU shares use its separately snapshotted ISR interval of 1999657 us.

```
isr_wake       1151.9/f 0.268/1.695/22.660 us 3.12%
  isr_pack       144.0/f 7.148/8.453/13.583 us 1.94%
  isr_dma_submit 144.0/f 0.576/0.923/6.965 us 0.21%
```

- Mean pack 8.453 us/call versus submit 0.923 us/call (9.2×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.12%, corresponding to 60.550 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.393×; Peak live show render: 0.746×; Show end: 0.448×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `grd_render` — 74.2% of the frame, 46.321 ms/f.
2. `grd_shader_draw` — 43.2% of the frame, 26.978 ms/f.
3. `grd_physics` — 8.3% of the frame, 5.169 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the effect-local physics/shading functions and shared scan/filter kernels. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
