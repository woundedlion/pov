# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_mindsplatter_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh MindSplatter profile_o3 110 16`). Raw capture:
`build/prof/mindsplatter_o3.log`, captured 2026-08-26 01:49 local. This replaces `profile_mindsplatter_teensy_2026-08-07.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile_o3 110 16` |

Image size: `FLASH: code:88,416, data:545,496, headers:9,160` / `RAM1: variables:315,296, code:62,648, padding:2,888, free:143,456` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 1473–1488
root counter cycles ÷ 600 MHz match the measured wall sum within **2.6 ppm**
(`tools/parse_profile.py build/prof/mindsplatter_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `msp_draw_particles` avg
22.60 ms/f, worst window 53.40 ms/f
(frames 1473–1488),
peak frame render 66.86 ms, spilled 10/1712
frames (0.6%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
The pass is yellow: 0.6% of frames exceed one display interval. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 1473–1488)

```text
frame                         101.49ms  60.89Mcyc 100.0%
  msp_draw_particles           53.40ms  32.04Mcyc  52.6%
    msp_particle_scan          53.40ms  32.04Mcyc  52.6%
      plot_ps_raster           39.14ms  23.48Mcyc  38.6% x915.1 43us/c
      plot_ps_deferred          3.17ms   1.90Mcyc   3.1% x915.1 3us/c
      plot_ps_gate              8.07ms   4.84Mcyc   8.0%
        plot_ps_cartesian_…    937.4us  562.4kcyc   0.9% x1367.6 1us/c
      plot_ps_tween             2.61ms   1.57Mcyc   2.6% x1367.6 2us/c
  msp_particle_step             4.94ms   2.96Mcyc   4.9% x1.0 4936us/c
  msp_timeline_step             49.6us   29.8kcyc   0.0% x1.0 50us/c
  canvas_clear                  86.6us   52.0kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           43.01ms  25.81Mcyc  42.4% x1.0 43013us/c
```

Wall min/avg/max = 24.23/101.49/126.26 ms. The `msp_draw_particles` scope averages 53.40 ms/f in this window, while the exact frame-render peak is 66.86 ms. That is 4.36 ms of overrun against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 1073–1088)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  msp_draw_particles            7.39ms   4.44Mcyc  11.8%
    msp_particle_scan           7.39ms   4.43Mcyc  11.8%
      plot_ps_raster            4.25ms   2.55Mcyc   6.8% x60.5 70us/c
      plot_ps_deferred         264.5us  158.7kcyc   0.4% x60.5 4us/c
      plot_ps_gate             983.5us  590.1kcyc   1.6%
        plot_ps_cartesian_…    462.7us  277.6kcyc   0.7% x880.8 1us/c
      plot_ps_tween             1.67ms   1.00Mcyc   2.7% x880.8 2us/c
  msp_particle_step             4.80ms   2.88Mcyc   7.7% x1.0 4801us/c
  msp_timeline_step             40.4us   24.3kcyc   0.1% x1.0 40us/c
  canvas_clear                  84.9us   50.9kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           50.10ms  30.06Mcyc  80.3% x1.0 50104us/c
```

Wall min/avg/max = 59.94/62.42/64.71 ms. The `msp_draw_particles` scope averages 7.39 ms/f in this window, while the exact frame-render peak is 13.77 ms. That is 48.73 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `msp_draw_particles` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `2` | — | 20/20 | — | 53.40 | 66.86 | 8.0 |
| `3` | — | 17/17 | — | 44.80 | 58.78 | 16.0 |
| `7` | — | 10/10 | — | 44.44 | 56.67 | 16.0 |
| `6` | — | 10/10 | — | 35.16 | 45.95 | 16.0 |
| `4` | — | 10/10 | — | 30.16 | 35.80 | 16.0 |
| `1` | — | 10/10 | — | 25.07 | 29.74 | 16.0 |
| `8` | — | 10/10 | — | 24.55 | 30.80 | 16.0 |
| `0` | — | 10/10 | — | 17.44 | 20.62 | 16.0 |
| `5` | — | 10/10 | — | 15.84 | 20.50 | 16.0 |

97 advance markers visit 8/8 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1871.7/frame  min/avg/max 0.36/1.54/12.20 us  cpu 2.84%
isr_pack          233.9/frame  min/avg/max 5.99/6.87/9.29 us  cpu 1.58%
isr_dma_submit    233.9/frame  min/avg/max 0.61/0.93/4.89 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.63% CPU. The worst render requires 7.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `msp_draw_particles` — 36.0% of measured root time, 22.60 ms/f average.
2. `msp_particle_scan` — 36.0% of measured root time, 22.60 ms/f average.
3. `plot_ps_raster` — 25.0% of measured root time, 15.71 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MindSplatter profile_o3 110 16` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 67.05 ms versus
66.86 ms here: global O3 lowers the peak by 0.19 ms (0.3%). O3-image minus shipping-image
size deltas are **FLASH code +22,760 B** and **ITCM code
+20,368 B**. Spill fractions are compared rather than raw counts:
shipping 39/1696 (2.3%)
versus O3 10/1712 (0.6%).
