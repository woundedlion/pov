# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_hyperlattice_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh HyperLattice profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"`). Raw capture:
`build/prof/hyperlattice_o3.log`, captured 2026-08-26 02:37 local. This replaces `profile_hyperlattice_teensy_2026-08-25.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, flags `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size: `FLASH: code:72,384, data:147,620, headers:8,348` / `RAM1: variables:315,040, code:30,040, padding:2,728, free:176,480` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 849–864
root counter cycles ÷ 600 MHz match the measured wall sum within **2.6 ppm**
(`tools/parse_profile.py build/prof/hyperlattice_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hl_shader_draw` avg
37.05 ms/f, worst window 44.97 ms/f
(frames 849–864),
peak frame render 49.92 ms, spilled 0/2688
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 545–560)

```text
frame                          62.08ms  37.25Mcyc 100.0%
  hl_shader_draw               43.16ms  25.90Mcyc  69.5% x1.0 43165us/c
  hl_timeline_step              11.8us    7.1kcyc   0.0% x1.0 12us/c
  canvas_clear                  87.1us   52.3kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           16.47ms   9.88Mcyc  26.5% x1.0 16471us/c
```

Wall min/avg/max = 55.95/62.08/67.22 ms. The `hl_shader_draw` scope averages 43.16 ms/f in this window, while the exact frame-render peak is 49.92 ms. That is 12.58 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 1233–1248)

```text
frame                          62.41ms  37.45Mcyc 100.0%
  hl_shader_draw               27.80ms  16.68Mcyc  44.5% x1.0 27799us/c
  hl_timeline_step               1.9us    1.2kcyc   0.0% x1.0 2us/c
  canvas_clear                  86.9us   52.2kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           32.31ms  19.39Mcyc  51.8% x1.0 32312us/c
```

Wall min/avg/max = 60.88/62.41/63.82 ms. The `hl_shader_draw` scope averages 27.80 ms/f in this window, while the exact frame-render peak is 31.08 ms. That is 31.42 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `hl_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `2` | — | 70/70 | — | 44.97 | 49.92 | 16.0 |
| `3` | — | 43/43 | — | 42.98 | 47.67 | 16.0 |
| `0` | — | 20/20 | — | 37.35 | 41.45 | 16.0 |
| `1` | — | 35/35 | — | 37.31 | 41.01 | 16.0 |

148 advance markers visit 3/3 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1145.3/frame  min/avg/max 0.42/1.62/25.14 us  cpu 2.99%
isr_pack          143.2/frame  min/avg/max 6.00/6.60/8.88 us  cpu 1.52%
isr_dma_submit    143.2/frame  min/avg/max 0.65/0.94/5.08 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.72% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `hl_shader_draw` — 59.4% of measured root time, 37.05 ms/f average.
2. `canvas_clear` — 0.1% of measured root time, 0.09 ms/f average.
3. `hl_timeline_step` — 0.0% of measured root time, 0.01 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=1600` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `20ca3cb48892795a4575dd9d16d31a699f79df75`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh HyperLattice profile_o3 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 51.19 ms versus
49.92 ms here: global O3 lowers the peak by 1.26 ms (2.5%). O3-image minus shipping-image
size deltas are **FLASH code +10,152 B** and **ITCM code
+7,728 B**. Spill fractions are compared rather than raw counts:
shipping 0/2688 (0%)
versus O3 0/2688 (0%).
