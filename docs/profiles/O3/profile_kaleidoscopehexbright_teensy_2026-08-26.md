# KaleidoscopeHexBright on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopehexbright_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeHexBright profile_o3 150 32`). Raw capture:
`build/prof/kaleidoscopehexbright_o3.log`, captured 2026-08-26 03:14 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexBright 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexBright profile_o3 150 32` |

Image size: `FLASH: code:83,960, data:151,932, headers:8,844` / `RAM1: variables:315,104, code:37,848, padding:27,688, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 385–416
root counter cycles ÷ 600 MHz match the measured wall sum within **1.4 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopehexbright_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
23.90 ms/f, worst window 26.49 ms/f
(frames 385–416),
peak frame render 32.41 ms, spilled 0/2368
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 385–416)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  fx_shader_draw               26.49ms  15.89Mcyc  42.4% x1.0 26486us/c
  fx_prepare_frame              3.65ms   2.19Mcyc   5.8% x1.0 3650us/c
  fx_advance                    2.02ms   1.21Mcyc   3.2% x1.0 2025us/c
  fx_timeline_step              44.3us   26.6kcyc   0.1% x1.0 44us/c
  canvas_clear                  86.2us   51.7kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           30.12ms  18.07Mcyc  48.3% x1.0 30125us/c
```

Wall min/avg/max = 62.20/62.43/62.66 ms. The `fx_shader_draw` scope averages 26.49 ms/f in this window, while the exact frame-render peak is 32.41 ms. That is 30.09 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 641–672)

```text
frame                          62.41ms  37.44Mcyc 100.0%
  fx_shader_draw               21.24ms  12.74Mcyc  34.0% x1.0 21237us/c
  fx_prepare_frame              3.33ms   2.00Mcyc   5.3% x1.0 3332us/c
  fx_advance                    2.03ms   1.22Mcyc   3.3% x1.0 2029us/c
  fx_timeline_step              90.4us   54.3kcyc   0.1% x1.0 90us/c
  canvas_clear                  85.7us   51.4kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           35.62ms  21.37Mcyc  57.1% x1.0 35620us/c
```

Wall min/avg/max = 61.46/62.41/63.06 ms. The `fx_shader_draw` scope averages 21.24 ms/f in this window, while the exact frame-render peak is 27.11 ms. That is 35.39 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `0` | — | 19/19 | — | 26.49 | 32.41 | 16.0 |
| `1` | — | 21/21 | — | 25.99 | 32.35 | 16.0 |
| `2` | — | 34/34 | — | 25.38 | 32.30 | 16.0 |

55 advance markers visit 2/2 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame  min/avg/max 0.39/1.57/24.45 us  cpu 2.90%
isr_pack          144.0/frame  min/avg/max 6.21/7.17/9.54 us  cpu 1.65%
isr_dma_submit    144.0/frame  min/avg/max 0.61/0.93/4.73 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.76% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 38.3% of measured root time, 23.90 ms/f average.
2. `fx_prepare_frame` — 6.0% of measured root time, 3.73 ms/f average.
3. `fx_advance` — 3.2% of measured root time, 2.03 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexBright`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeHexBright profile_o3 150 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 35.55 ms versus
32.41 ms here: global O3 lowers the peak by 3.14 ms (8.8%). O3-image minus shipping-image
size deltas are **FLASH code +14,520 B** and **ITCM code
+11,712 B**. Spill fractions are compared rather than raw counts:
shipping 0/2368 (0%)
versus O3 0/2368 (0%).
