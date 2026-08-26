# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_mobiusgrid_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh MobiusGrid profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"`). Raw capture:
`build/prof/mobiusgrid_o3.log`, captured 2026-08-26 01:30 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, flags `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh MobiusGrid profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size: `FLASH: code:83,456, data:152,404, headers:8,876` / `RAM1: variables:315,104, code:36,792, padding:28,744, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 1841–1856
root counter cycles ÷ 600 MHz match the measured wall sum within **2.0 ppm**
(`tools/parse_profile.py build/prof/mobiusgrid_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
18.00 ms/f, worst window 18.41 ms/f
(frames 1841–1856),
peak frame render 21.44 ms, spilled 0/2688
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 593–608)

```text
frame                          62.46ms  37.47Mcyc 100.0%
  fx_shader_draw               18.09ms  10.85Mcyc  29.0% x1.0 18089us/c
  fx_prepare_frame             434.7us  260.8kcyc   0.7% x1.0 435us/c
  fx_advance                    2.01ms   1.21Mcyc   3.2% x1.0 2012us/c
  fx_timeline_step              88.6us   53.1kcyc   0.1% x1.0 89us/c
  canvas_clear                  85.3us   51.2kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           41.74ms  25.04Mcyc  66.8% x1.0 41736us/c
```

Wall min/avg/max = 59.65/62.46/65.00 ms. The `fx_shader_draw` scope averages 18.09 ms/f in this window, while the exact frame-render peak is 21.44 ms. That is 41.06 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 2369–2384)

```text
frame                          62.41ms  37.45Mcyc 100.0%
  fx_shader_draw               17.78ms  10.67Mcyc  28.5% x1.0 17775us/c
  fx_prepare_frame             383.8us  230.3kcyc   0.6% x1.0 384us/c
  fx_advance                    2.22ms   1.33Mcyc   3.6% x1.0 2216us/c
  fx_timeline_step              52.6us   31.6kcyc   0.1% x1.0 53us/c
  canvas_clear                  85.2us   51.1kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           41.89ms  25.13Mcyc  67.1% x1.0 41885us/c
```

Wall min/avg/max = 62.20/62.41/62.51 ms. The `fx_shader_draw` scope averages 17.78 ms/f in this window, while the exact frame-render peak is 20.58 ms. That is 41.92 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `1` | — | 63/63 | — | 18.41 | 21.39 | 16.0 |
| `2` | — | 67/67 | — | 18.41 | 21.44 | 16.0 |
| `0` | — | 38/38 | — | 18.09 | 20.76 | 16.0 |

130 advance markers visit 2/2 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.9/frame  min/avg/max 0.40/1.50/14.78 us  cpu 2.76%
isr_pack          144.1/frame  min/avg/max 5.99/6.55/8.78 us  cpu 1.51%
isr_dma_submit    144.1/frame  min/avg/max 0.63/0.93/7.89 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.48% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 28.9% of measured root time, 18.00 ms/f average.
2. `fx_advance` — 3.5% of measured root time, 2.19 ms/f average.
3. `fx_prepare_frame` — 0.7% of measured root time, 0.42 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MobiusGrid`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MobiusGrid profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 22.37 ms versus
21.44 ms here: global O3 lowers the peak by 0.93 ms (4.2%). O3-image minus shipping-image
size deltas are **FLASH code +13,448 B** and **ITCM code
+10,560 B**. Spill fractions are compared rather than raw counts:
shipping 0/2688 (0%)
versus O3 0/2688 (0%).
