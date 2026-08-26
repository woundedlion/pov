# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopehexoil_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeHexOil profile_o3 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"`). Raw capture:
`build/prof/kaleidoscopehexoil_o3.log`, captured 2026-08-26 02:49 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 140 s capture, flags `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexOil profile_o3 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size: `FLASH: code:83,664, data:155,868, headers:8,276` / `RAM1: variables:315,104, code:36,584, padding:28,952, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 481–496
root counter cycles ÷ 600 MHz match the measured wall sum within **3.7 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopehexoil_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
32.64 ms/f, worst window 35.28 ms/f
(frames 481–496),
peak frame render 38.94 ms, spilled 0/2208
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 497–512)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  fx_shader_draw               35.18ms  21.11Mcyc  56.4% x1.0 35184us/c
  fx_prepare_frame             889.8us  533.9kcyc   1.4% x1.0 890us/c
  fx_advance                    2.27ms   1.36Mcyc   3.6% x1.0 2268us/c
  fx_timeline_step              50.4us   30.2kcyc   0.1% x1.0 50us/c
  canvas_clear                  87.5us   52.5kcyc   0.1% x1.0 88us/c
  canvas_buffer_wait           23.93ms  14.36Mcyc  38.3% x1.0 23932us/c
```

Wall min/avg/max = 61.43/62.43/63.41 ms. The `fx_shader_draw` scope averages 35.18 ms/f in this window, while the exact frame-render peak is 38.94 ms. That is 23.56 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 641–656)

```text
frame                          62.39ms  37.43Mcyc 100.0%
  fx_shader_draw               30.43ms  18.26Mcyc  48.8% x1.0 30430us/c
  fx_prepare_frame             851.2us  510.8kcyc   1.4% x1.0 851us/c
  fx_advance                    2.29ms   1.37Mcyc   3.7% x1.0 2286us/c
  fx_timeline_step              95.9us   57.6kcyc   0.2% x1.0 96us/c
  canvas_clear                  86.9us   52.2kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           28.62ms  17.17Mcyc  45.9% x1.0 28622us/c
```

Wall min/avg/max = 61.48/62.39/63.15 ms. The `fx_shader_draw` scope averages 30.43 ms/f in this window, while the exact frame-render peak is 34.17 ms. That is 28.33 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `0` | — | 38/38 | — | 35.28 | 38.94 | 16.0 |
| `2` | — | 67/67 | — | 34.23 | 38.41 | 16.0 |
| `1` | — | 33/33 | — | 34.04 | 38.65 | 16.0 |

100 advance markers visit 2/2 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.3/frame  min/avg/max 0.39/1.61/17.84 us  cpu 2.96%
isr_pack          144.0/frame  min/avg/max 6.20/7.28/9.33 us  cpu 1.67%
isr_dma_submit    144.0/frame  min/avg/max 0.60/0.93/4.85 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.84% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 52.3% of measured root time, 32.64 ms/f average.
2. `fx_advance` — 3.6% of measured root time, 2.23 ms/f average.
3. `fx_prepare_frame` — 1.1% of measured root time, 0.71 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=1600` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexOil`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeHexOil profile_o3 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 39.04 ms versus
38.94 ms here: global O3 lowers the peak by 0.10 ms (0.3%). O3-image minus shipping-image
size deltas are **FLASH code +13,456 B** and **ITCM code
+10,608 B**. Spill fractions are compared rather than raw counts:
shipping 0/2208 (0%)
versus O3 0/2208 (0%).
