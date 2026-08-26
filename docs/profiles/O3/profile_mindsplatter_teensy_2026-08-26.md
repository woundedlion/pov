# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_mindsplatter_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh MindSplatter profile_o3 110 16`). Raw capture:
`build/prof/mindsplatter_o3.log`, captured 2026-08-26 07:45 local. This refresh replaces the
earlier 01:49 capture with the final clean cadence-reclaim image.
Its `.provenance`, `_envdump.txt`, and `_build.log` sidecars identify the exact
source, compiler flags, ELF hashes, and image sizes used below.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `0df961b818ae08c5f58639edb93d812164a9356a` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile_o3 110 16` |

Image size: `FLASH: code:88,416, data:545,504, headers:9,152` / `RAM1: variables:315,296, code:62,648, padding:2,888, free:143,456` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 993–1008 root counter
cycles ÷ 600 MHz match the measured wall sum within **0.3 ppm**
(`tools/parse_profile.py build/prof/mindsplatter_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `msp_draw_particles` avg
19.22 ms/f, worst window 42.44 ms/f
(frames 993–1008), peak frame render 52.79 ms,
spilled 0/1728 frames (0.0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
The pass is green: every captured frame stays within one display interval, and
the worst frame retains 9.71 ms of render margin. The
`canvas_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: marker-defined presets own their following transitions; the
first section is the worst exact-render window and the second is the
lowest-cost marked regime.

### Worst exact-render regime (window frames 993–1008)

```text
frame                         62.24ms 37.35Mcyc 100.0%
  msp_draw_particles          42.44ms 25.46Mcyc  68.2%
    msp_particle_scan         42.43ms 25.46Mcyc  68.2%
      plot_ps_raster          32.59ms 19.55Mcyc  52.4% x540.7 60us/c
      plot_ps_deferred        435.2us 261.1kcyc   0.7% x540.7 483cyc/c
      plot_ps_gate             5.82ms  3.49Mcyc   9.3%
        plot_ps_cartesian_g…   1.02ms 612.6kcyc   1.6% x1536.9 399cyc/c
      plot_ps_tween            3.18ms  1.91Mcyc   5.1% x1536.9 2us/c
  msp_particle_step            2.88ms  1.73Mcyc   4.6% x1.0 2883us/c
  msp_timeline_step            43.4us  26.0kcyc   0.1% x1.0 43us/c
  canvas_clear                 88.7us  53.2kcyc   0.1% x1.0 89us/c
  canvas_buffer_wait          16.79ms 10.07Mcyc  27.0% x1.0 16790us/c
```

Wall min/avg/max = 50.70/62.24/73.43 ms.
`msp_draw_particles` averages 42.44 ms/f here, while the exact frame-render
peak is 52.79 ms. The phase retains 9.71 ms against
one display interval; display-sync wait is idle rather than render work.

### Lowest-cost marked regime (window frames 257–272)

```text
frame                         62.08ms 37.25Mcyc 100.0%
  msp_draw_particles           6.08ms  3.65Mcyc   9.8%
    msp_particle_scan          6.08ms  3.65Mcyc   9.8%
      plot_ps_raster           2.60ms  1.56Mcyc   4.2% x52.2 50us/c
      plot_ps_deferred         54.6us  32.7kcyc   0.1% x52.2 1us/c
      plot_ps_gate             1.09ms 655.9kcyc   1.8%
        plot_ps_cartesian_g…  567.7us 340.6kcyc   0.9% x1043.0 327cyc/c
      plot_ps_tween            2.09ms  1.25Mcyc   3.4% x1043.0 2us/c
  msp_particle_step            2.07ms  1.24Mcyc   3.3% x1.0 2068us/c
  msp_timeline_step            39.8us  23.9kcyc   0.1% x1.0 40us/c
  canvas_clear                 85.3us  51.2kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait          53.80ms 32.28Mcyc  86.7% x1.0 53803us/c
```

Wall min/avg/max = 56.72/62.08/68.00 ms.
`msp_draw_particles` falls to 6.08 ms/f and the exact render peak is
12.65 ms, leaving 49.85 ms of margin. The spread isolates
the preset-dependent particle raster and gate work from the fixed display cadence.

### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `msp_draw_particles` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `7` | — | 10/10 | — | 42.44 | 52.79 | 16.0 |
| `2` | — | 20/20 | — | 35.52 | 47.58 | 16.0 |
| `6` | — | 10/10 | — | 35.18 | 45.94 | 16.0 |
| `4` | — | 10/10 | — | 25.05 | 29.09 | 16.0 |
| `1` | — | 10/10 | — | 24.36 | 28.29 | 16.0 |
| `8` | — | 10/10 | — | 23.25 | 29.02 | 16.0 |
| `3` | — | 18/18 | — | 23.07 | 29.14 | 16.0 |
| `0` | — | 10/10 | — | 17.45 | 20.64 | 16.0 |
| `5` | — | 10/10 | — | 15.82 | 19.68 | 16.0 |

All 8 presets are visited and the sequence wraps back to preset 0.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1148.1/frame  min/avg/max 0.37/1.59/11.89 us  cpu 2.93%
isr_pack         143.5/frame  min/avg/max 5.99/6.98/9.32 us  cpu 1.60%
isr_dma_submit   143.5/frame  min/avg/max 0.63/0.93/1.01 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are
  reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains
  isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.74% CPU. The worst render
  already fits one display interval with 9.71 ms of margin.

## Summary ranking

1. `msp_draw_particles` — 30.8% of measured root time, 19.22 ms/f average.
2. `msp_particle_scan` — 30.8% of measured root time, 19.22 ms/f average.
3. `plot_ps_raster` — 21.8% of measured root time, 13.60 ms/f average.

No matched WASM/native datum is asserted here; the final paired on-device
shipping capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under whichever scope first enters it; its subtree is
  hidden in windows where that parent has zero calls, and calls approximate
  blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths
  identified by the counter tree; global O3 compiles the entire single-effect
  image.
- No dwell-compression or ordered-cycle override was used.
- Presets 2/8, 3/8, and 7/8 select Cube geometry. Emitters remain uncapped and
  every active emitter continues to spawn each frame.
- Provenance records clean O3 source `0df961b818ae08c5f58639edb93d812164a9356a` and clean paired
  shipping source `0df961b818ae08c5f58639edb93d812164a9356a`; the capture artifacts record no
  uncommitted source state.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MindSplatter`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MindSplatter profile_o3 110 16` builds, flashes, and captures the
global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The final paired shipping capture peaks at 52.77 ms versus
52.79 ms here: global O3 raises the peak by 0.02 ms
(0.04%). O3-image minus shipping-image size deltas are **FLASH code
+22,800 B** and **ITCM code +20,400 B**. Both captures are green:
shipping 0/1728 spills versus O3
0/1728.
