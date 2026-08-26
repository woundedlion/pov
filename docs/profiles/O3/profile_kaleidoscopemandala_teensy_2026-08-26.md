# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopemandala_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeMandala profile_o3 150 32`). Raw capture:
`build/prof/kaleidoscopemandala_o3.log`, captured 2026-08-26 03:17 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeMandala 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeMandala profile_o3 150 32` |

Image size: `FLASH: code:83,528, data:152,080, headers:9,128` / `RAM1: variables:315,104, code:37,608, padding:27,928, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 1153–1184
root counter cycles ÷ 600 MHz match the measured wall sum within **4.0 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopemandala_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
28.17 ms/f, worst window 30.85 ms/f
(frames 1153–1184),
peak frame render 36.86 ms, spilled 0/2368
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 1921–1952)

```text
frame                          62.41ms  37.45Mcyc 100.0%
  fx_shader_draw               30.11ms  18.07Mcyc  48.2% x1.0 30112us/c
  fx_prepare_frame              3.90ms   2.34Mcyc   6.2% x1.0 3899us/c
  fx_advance                    2.21ms   1.33Mcyc   3.5% x1.0 2215us/c
  fx_timeline_step              97.2us   58.4kcyc   0.2% x1.0 97us/c
  canvas_clear                  86.6us   52.0kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           25.99ms  15.59Mcyc  41.6% x1.0 25990us/c
```

Wall min/avg/max = 61.03/62.41/63.63 ms. The `fx_shader_draw` scope averages 30.11 ms/f in this window, while the exact frame-render peak is 36.86 ms. That is 25.64 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 705–736)

```text
frame                          62.40ms  37.44Mcyc 100.0%
  fx_shader_draw               26.34ms  15.80Mcyc  42.2% x1.0 26337us/c
  fx_prepare_frame              3.91ms   2.34Mcyc   6.3% x1.0 3905us/c
  fx_advance                    2.27ms   1.36Mcyc   3.6% x1.0 2266us/c
  fx_timeline_step              99.7us   59.8kcyc   0.2% x1.0 100us/c
  canvas_clear                  86.4us   51.9kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           29.69ms  17.81Mcyc  47.6% x1.0 29689us/c
```

Wall min/avg/max = 61.53/62.40/63.22 ms. The `fx_shader_draw` scope averages 26.34 ms/f in this window, while the exact frame-render peak is 33.28 ms. That is 29.22 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `2` | — | 34/34 | — | 30.85 | 36.79 | 16.0 |
| `1` | — | 21/21 | — | 30.11 | 36.86 | 16.0 |
| `0` | — | 19/19 | — | 29.61 | 33.57 | 16.0 |

55 advance markers visit 2/2 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.9/frame  min/avg/max 0.40/1.59/25.41 us  cpu 2.94%
isr_pack          144.0/frame  min/avg/max 6.29/7.07/9.39 us  cpu 1.62%
isr_dma_submit    144.0/frame  min/avg/max 0.66/0.94/11.31 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.77% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 45.1% of measured root time, 28.17 ms/f average.
2. `fx_advance` — 3.6% of measured root time, 2.24 ms/f average.
3. `fx_prepare_frame` — 3.3% of measured root time, 2.05 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `20ca3cb48892795a4575dd9d16d31a699f79df75` and paired shipping source
  `20ca3cb48892795a4575dd9d16d31a699f79df75`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeMandala`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeMandala profile_o3 150 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 40.29 ms versus
36.86 ms here: global O3 lowers the peak by 3.43 ms (8.5%). O3-image minus shipping-image
size deltas are **FLASH code +13,656 B** and **ITCM code
+11,440 B**. Spill fractions are compared rather than raw counts:
shipping 0/2368 (0%)
versus O3 0/2368 (0%).
