# LatticeMelt on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_latticemelt_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh LatticeMelt profile_o3 230 16 "-D HS_PROFILE_EPOCH_REVS=2400"`). Raw capture:
`build/prof/latticemelt_o3.log`, captured 2026-08-26 03:22 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | LatticeMelt 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 230 s capture, flags `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh LatticeMelt profile_o3 230 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size: `FLASH: code:87,448, data:151,948, headers:8,412` / `RAM1: variables:315,104, code:38,088, padding:27,448, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 369–384
root counter cycles ÷ 600 MHz match the measured wall sum within **1.3 ppm**
(`tools/parse_profile.py build/prof/latticemelt_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
93.32 ms/f, worst window 99.82 ms/f
(frames 369–384),
peak frame render 104.75 ms, spilled 1824/1824
frames (100%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
The pass is red: 100% of frames exceed one display interval. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 385–400)

```text
frame                         124.78ms  74.87Mcyc 100.0%
  fx_shader_draw               99.11ms  59.47Mcyc  79.4% x1.0 99115us/c
  fx_prepare_frame              1.06ms  637.1kcyc   0.9% x1.0 1062us/c
  fx_advance                    2.22ms   1.33Mcyc   1.8% x1.0 2221us/c
  fx_timeline_step              50.3us   30.2kcyc   0.0% x1.0 50us/c
  canvas_clear                  86.3us   51.8kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           22.23ms  13.34Mcyc  17.8% x1.0 22228us/c
```

Wall min/avg/max = 122.72/124.78/126.70 ms. The `fx_shader_draw` scope averages 99.11 ms/f in this window, while the exact frame-render peak is 104.75 ms. That is 42.25 ms of overrun against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 785–800)

```text
frame                         125.06ms  75.04Mcyc 100.0%
  fx_shader_draw               86.18ms  51.71Mcyc  68.9% x1.0 86176us/c
  fx_prepare_frame             931.0us  558.6kcyc   0.7% x1.0 931us/c
  fx_advance                    2.25ms   1.35Mcyc   1.8% x1.0 2246us/c
  fx_timeline_step             111.1us   66.7kcyc   0.1% x1.0 111us/c
  canvas_clear                  85.5us   51.3kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           35.50ms  21.30Mcyc  28.4% x1.0 35495us/c
```

Wall min/avg/max = 122.28/125.06/127.08 ms. The `fx_shader_draw` scope averages 86.18 ms/f in this window, while the exact frame-render peak is 91.09 ms. That is 28.59 ms of overrun against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `0` | — | 38/38 | — | 99.82 | 104.75 | 8.0 |
| `2` | — | 67/67 | — | 99.11 | 103.80 | 8.0 |
| `1` | — | 9/9 | — | 87.90 | 98.07 | 8.0 |

76 advance markers visit 2/2 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         2301.9/frame  min/avg/max 0.39/9.95/36.28 us  cpu 18.35%
isr_pack          287.8/frame  min/avg/max 5.99/6.81/9.08 us  cpu 1.56%
isr_dma_submit    287.8/frame  min/avg/max 0.73/0.93/9.58 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 20.12% CPU. The worst render requires 67.6% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 74.7% of measured root time, 93.32 ms/f average.
2. `fx_advance` — 1.8% of measured root time, 2.23 ms/f average.
3. `fx_prepare_frame` — 0.7% of measured root time, 0.84 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=2400` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `20ca3cb48892795a4575dd9d16d31a699f79df75` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=LatticeMelt`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh LatticeMelt profile_o3 230 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 43.73 ms versus
104.75 ms here: global O3 raises the peak by 61.02 ms (139.6%). O3-image minus shipping-image
size deltas are **FLASH code +16,592 B** and **ITCM code
+11,952 B**. Spill fractions are compared rather than raw counts:
shipping 0/1728 (0%)
versus O3 1824/1824 (100%).
