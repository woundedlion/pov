# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-11-17, **selective -O3**)

Point-in-time snapshot of the authored show slot. Reproduce with the scoped command in Setup.
Raw capture: `C:\work\temp\holosphere-profile-results-20261117\build\prof\displacementfield_ship.log`. Replaces `profile_displacementfield_teensy_2026-11-12.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; live flywheel and DMA ISRs |
| Image | `profile`; shipping `-Os` base, newlib-nano; selective `HS_O3`: the effect-local ball/noise field evaluators and `core/render/scan/ring_stack.h` |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, source `6ac5a4f049fd06562df1530888d9a50fcbfa9f21` |
| Method | Show scope = 121 s / 1936 nominal frames at 16 fps; window 16. Raw capture 190 s; later frames excluded. Profiling epoch 3600 s; no reset in the capture. |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh DisplacementField profile 124 16`; analyze only the 121 s show budget (clean checkout of the source above) |

Profile image size: `FLASH: code:118168, data:150728, headers:8608; free:1754112 / RAM1: variables:314464, code:25320, padding:7448; free:177056 / RAM2: variables:520064; free:4224`.

Shipping Phantasm image size: `FLASH: code:1167656, data:430424, headers:8576; free:424960 / RAM1: variables:313760, code:172072, padding:24536; free:13920 / RAM2: variables:520064; free:4224`.

Exactness cross-check: frames 593–608 root cycles / 600 MHz = 999583.785 us; measured wall sum 999584 us, delta 0.22 ppm. `parse_profile.py ... validate` passes frame coverage, configuration and no-reset checks.

Captured 2026-11-17 08:38, America/Los_Angeles. Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256: `164c33340a3d655ee8c01588aaef2235e03702787ad2596c6a971f47ee5b39ca`. The metrics sidecar retains verified ELF/envdump hashes and immutable artifact paths.

## Frame cadence

**Actual-show scope:** the playlist assigns 121 s (1936 nominal frames at 16 fps). Show-frame peak render is **44.616 ms**, spilled **0/1936**. Setup frame 1 is **18.450 ms** and remains separate. Live show peak is **44.616 ms**, spilled **0/1935**.

At 480 RPM a half revolution is 62.5 ms, giving 16 fps. This effect renders the clipped quadrant, containing 10,368 pixels. Render telemetry subtracts `canvas_buffer_wait`, the round-up idle until the next display flip; wall maxima are not render peaks.

The ranked index uses only show frames 2–1936. The raw log retains 3016 frame records; the 1080 records after the show budget are excluded from peaks, spill denominators, phase readouts and per-preset costs. All transitions inside the show budget remain included. Counter trees use complete in-scope windows.

README cells: peak 🟢 44.616 (2), spilled 🟢 0/1935 (0.00%).

README setup frame 1: 18.450 ms.

## Phase-by-phase readout

Representative windows distinguish setup, the largest individual live-show render, qualified preset cost when applicable, and the end of the show slot. Every window is inside the show budget. Tree entries are inclusive per-frame cycle measurements; window means do not replace individual-frame peaks.

### Startup (frames 1–16)

```
frame                    57.776 ms  34.666 Mcyc 100.0%
  pov_preserve_half      135.30 us   81.18 kcyc   0.2% x0.9375 144 us/call
  df_timeline_step       25.699 ms  15.419 Mcyc  44.5%
    df_draw_rings        25.592 ms  15.355 Mcyc  44.3%
      df_fused_scan       9.973 ms   5.984 Mcyc  17.3% x0.5625 17729 us/call
      df_lut_bake         6.822 ms   4.093 Mcyc  11.8%
        df_hue_table_prep  199.23 us  119.54 kcyc   0.3% x8.0625 25 us/call
        df_octave_noise   4.695 ms   2.817 Mcyc   8.1% x26.8125 175 us/call
      df_chunk_cull      284.35 us  170.61 kcyc   0.5% x26.875 11 us/call
      df_lut_bake         1.546 ms  927.77 kcyc   2.7%
        df_hue_table_prep   35.54 us   21.33 kcyc   0.1% x5.25 7 us/call
      df_chunk_cull      349.09 us  209.45 kcyc   0.6% x18 19 us/call
      df_fused_scan       4.886 ms   2.932 Mcyc   8.5%
        ring_stack_table   1.617 ms  970.39 kcyc   2.8% x1.5 1078 us/call
        ring_stack_table   86.44 us   51.87 kcyc   0.1% x0.0625 1383 us/call
    animation_noise_clock    0.17 us    0.10 kcyc   0.0% x1 0 us/call
  canvas_clear            96.71 us   58.03 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait     31.823 ms  19.094 Mcyc  55.1% x1 31823 us/call
  df_prepare_fields        1.86 us    1.11 kcyc   0.0% x1 2 us/call
```

Wall min/avg/max = 17.602/57.776/77.507 ms. `df_timeline_step` averages 25.699 ms/f; total render averages 25.953 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Peak live show render (frames 593–608)

```
ring_stack_table          1.903 ms   1.142 Mcyc   3.0% x2 951 us/call
frame                    62.474 ms  37.484 Mcyc 100.0%
  pov_preserve_half      142.16 us   85.29 kcyc   0.2% x1 142 us/call
  df_timeline_step       42.991 ms  25.795 Mcyc  68.8%
    df_draw_rings        42.849 ms  25.709 Mcyc  68.6%
      df_fused_scan      25.638 ms  15.383 Mcyc  41.0% x1 25638 us/call
      df_lut_bake        15.029 ms   9.017 Mcyc  24.1%
        df_hue_table_prep   1.134 ms  680.52 kcyc   1.8% x21.625 52 us/call
        df_octave_noise  10.377 ms   6.226 Mcyc  16.6% x44.1875 235 us/call
      df_chunk_cull      506.26 us  303.75 kcyc   0.8% x45.875 11 us/call
    animation_noise_clock    0.17 us    0.10 kcyc   0.0% x1 0 us/call
  canvas_clear            88.77 us   53.26 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     19.215 ms  11.529 Mcyc  30.8% x1 19215 us/call
  df_prepare_fields        1.49 us    0.90 kcyc   0.0% x1 1 us/call
```

Wall min/avg/max = 60.646/62.474/64.331 ms. `df_timeline_step` averages 42.991 ms/f; total render averages 43.259 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Transition / non-hold regime (frames 737–752)

```
ring_stack_table          1.897 ms   1.138 Mcyc   3.0% x2 948 us/call
frame                    62.384 ms  37.430 Mcyc 100.0%
  pov_preserve_half      141.93 us   85.16 kcyc   0.2% x1 142 us/call
  df_timeline_step       42.772 ms  25.663 Mcyc  68.6%
    df_draw_rings        42.631 ms  25.579 Mcyc  68.3%
      df_fused_scan      25.688 ms  15.413 Mcyc  41.2% x1 25688 us/call
      df_lut_bake        14.785 ms   8.871 Mcyc  23.7%
        df_hue_table_prep   1.083 ms  649.65 kcyc   1.7% x22.1875 49 us/call
        df_octave_noise  10.277 ms   6.166 Mcyc  16.5% x43.5 236 us/call
      df_chunk_cull      495.49 us  297.30 kcyc   0.8% x45.5 11 us/call
    animation_noise_clock    0.17 us    0.10 kcyc   0.0% x1 0 us/call
  canvas_clear            89.37 us   53.62 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     19.345 ms  11.607 Mcyc  31.0% x1 19345 us/call
  df_prepare_fields        1.46 us    0.88 kcyc   0.0% x1 1 us/call
```

Wall min/avg/max = 61.175/62.383/64.033 ms. `df_timeline_step` averages 42.772 ms/f; total render averages 43.039 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Show end (frames 1921–1936)

```
ring_stack_table          1.227 ms  736.06 kcyc   2.0% x2 613 us/call
frame                    62.490 ms  37.494 Mcyc 100.0%
  pov_preserve_half      145.91 us   87.54 kcyc   0.2% x1 146 us/call
  df_timeline_step       29.307 ms  17.584 Mcyc  46.9%
    df_draw_rings        29.163 ms  17.498 Mcyc  46.7%
      df_fused_scan      16.963 ms  10.178 Mcyc  27.1% x1 16963 us/call
      df_lut_bake        10.392 ms   6.235 Mcyc  16.6%
        df_hue_table_prep  522.36 us  313.42 kcyc   0.8% x19.1875 27 us/call
        df_octave_noise   7.421 ms   4.453 Mcyc  11.9% x32.125 231 us/call
      df_chunk_cull      433.87 us  260.32 kcyc   0.7% x39.0625 11 us/call
    animation_noise_clock    0.17 us    0.10 kcyc   0.0% x1 0 us/call
  canvas_clear            88.81 us   53.29 kcyc   0.1% x1 89 us/call
  canvas_buffer_wait     32.912 ms  19.747 Mcyc  52.7% x1 32912 us/call
  df_prepare_fields        1.59 us    0.95 kcyc   0.0% x1 2 us/call
```

Wall min/avg/max = 61.637/62.489/63.311 ms. `df_timeline_step` averages 29.307 ms/f; total render averages 29.578 ms/f. Display synchronization accounts for the remaining buffer-wait time.

### Per-preset table

Noise has a 150-frame ramp-in, a 600-frame full-gain dwell and a 150-frame ramp-out. The noise row uses a full-gain window. The balls row uses a window from the 900-frame spawning phase, when the ball population changes. Cadence includes draining and phase handoffs.

The show-scoped data covers 2 preset identities. Completing the scheduled show is the coverage criterion; an additional wrap is not required. Cadence buckets include in-show transitions; qualified windows and live peaks are restricted to the same show scope. The fps column is the worst individual-frame cadence across each preset's qualified show windows.

| Preset | Qualified windows | Representative frames | Selection | V/E/F/I | Blended px/f | Draw ms/f | Render ms/f | fps | Live peak ms | Live spilled |
|---|--:|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 · noise | 36 | 593–608 | full-gain noise window | — | — | 42.991 | 43.259 | 16 | 44.616 | 0/923 |
| 2 · balls | 55 | 1617–1632 | ball-spawning regime window | — | — | 36.219 | 36.523 | 16 | 37.404 | 0/1012 |

### Per-pixel figures

This image emits no complete blended-pixel counter for the direct render path. Blended coverage and cycles per blend cannot be inferred from draw-call counts. The clipped quadrant has 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Representative peak window, frames 593–608; CPU shares use its separately snapshotted ISR interval of 1001844 us.

```
isr_wake       1154.2/f 0.285/1.661/28.775 us 3.06%
  isr_pack       144.2/f 6.845/7.993/10.773 us 1.84%
  isr_dma_submit 144.2/f 0.601/0.910/1.728 us 0.20%
```

- Mean pack 7.993 us/call versus submit 0.910 us/call (8.8×).
- Pack and DMA submit measure CPU work; submission does not measure the asynchronous wire transfer.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, framing gaps are included in the 230 us asynchronous transfer budget.
- Inclusive wake share is 3.06%, corresponding to 60.588 ms foreground capacity per display window. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured.

The inclusive render-time ratio to the 62.5 ms budget is Startup: 0.415×; Peak live show render: 0.692×; Transition / non-hold regime: 0.689×; Show end: 0.473×. Ratios above 1 require that speedup; ratios below 1 fit. These render measurements already include interrupts, so the wake share is not subtracted again.

## Summary ranking

1. `df_timeline_step` — 68.8% of the frame, 42.991 ms/f.
2. `df_draw_rings` — 68.6% of the frame, 42.849 ms/f.
3. `df_fused_scan` — 41.0% of the frame, 25.638 ms/f.

## Caveats

- Foreground scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parenting depends on its first caller; tagged MIXED-PARENT or DUPLICATE-NAME rows are non-exclusive and must not be summed with siblings.
- Existing instrumentation overhead remains included; no per-pixel scopes, PC sampling, replay, or count-only overrides were added for timing.
- Shipping selective-O3 regions active on this path: the effect-local ball/noise field evaluators and `core/render/scan/ring_stack.h`. This is the shipping `-Os` configuration.
- Full cycle means the effect's authored show slot. Authored preset speed and transitions are unchanged; capture padding and later repetitions are excluded from every reported show cost.
- This report applies the playlist's nominal frame budget to the single-effect capture. It does not directly measure the wall-clock handoff between effects in the full playlist.
- The measured source tree was clean. Source/compiler/ELF/envdump attestations identify the exact images; board provenance is recorded in Setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=16`. Use the supported locked wrapper with the explicit show-sized duration and window in Reproduce. Acquisition margin is excluded from the report.
