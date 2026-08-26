# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopeflowers_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeFlowers profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"`). Raw capture:
`build/prof/kaleidoscopeflowers_o3.log`, captured 2026-08-26 03:03 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, flags `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeFlowers profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size: `FLASH: code:84,752, data:151,984, headers:9,024` / `RAM1: variables:315,104, code:37,848, padding:27,688, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 561–576
root counter cycles ÷ 600 MHz match the measured wall sum within **1.6 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopeflowers_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
27.61 ms/f, worst window 29.63 ms/f
(frames 561–576),
peak frame render 33.38 ms, spilled 0/4128
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 1905–1920)

```text
frame                          62.42ms  37.45Mcyc 100.0%
  fx_shader_draw               29.11ms  17.47Mcyc  46.6% x1.0 29115us/c
  fx_prepare_frame              1.24ms  741.8kcyc   2.0% x1.0 1236us/c
  fx_advance                    2.01ms   1.21Mcyc   3.2% x1.0 2010us/c
  fx_timeline_step             102.8us   61.7kcyc   0.2% x1.0 103us/c
  canvas_clear                  85.8us   51.5kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           29.86ms  17.91Mcyc  47.8% x1.0 29857us/c
```

Wall min/avg/max = 60.75/62.42/63.91 ms. The `fx_shader_draw` scope averages 29.11 ms/f in this window, while the exact frame-render peak is 33.38 ms. That is 29.12 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 721–736)

```text
frame                          62.39ms  37.44Mcyc 100.0%
  fx_shader_draw               25.75ms  15.45Mcyc  41.3% x1.0 25747us/c
  fx_prepare_frame             412.2us  247.3kcyc   0.7% x1.0 412us/c
  fx_advance                    2.00ms   1.20Mcyc   3.2% x1.0 1997us/c
  fx_timeline_step              94.9us   57.0kcyc   0.2% x1.0 95us/c
  canvas_clear                  85.8us   51.5kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           34.04ms  20.43Mcyc  54.6% x1.0 34044us/c
```

Wall min/avg/max = 62.02/62.39/62.76 ms. The `fx_shader_draw` scope averages 25.75 ms/f in this window, while the exact frame-render peak is 28.52 ms. That is 33.98 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `0` | — | 38/38 | — | 29.63 | 32.67 | 16.0 |
| `2` | — | 85/85 | — | 29.60 | 32.83 | 16.0 |
| `3` | — | 68/68 | — | 29.57 | 33.38 | 16.0 |
| `1` | — | 67/67 | — | 29.55 | 33.11 | 16.0 |

220 advance markers visit 3/3 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame  min/avg/max 0.39/1.56/11.54 us  cpu 2.88%
isr_pack          144.0/frame  min/avg/max 6.16/7.15/9.61 us  cpu 1.64%
isr_dma_submit    144.0/frame  min/avg/max 0.72/0.93/1.51 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.73% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 44.3% of measured root time, 27.61 ms/f average.
2. `fx_advance` — 3.2% of measured root time, 1.99 ms/f average.
3. `fx_prepare_frame` — 1.3% of measured root time, 0.81 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=2400` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeFlowers`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeFlowers profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 36.69 ms versus
33.38 ms here: global O3 lowers the peak by 3.31 ms (9.0%). O3-image minus shipping-image
size deltas are **FLASH code +14,576 B** and **ITCM code
+11,712 B**. Spill fractions are compared rather than raw counts:
shipping 0/4128 (0%)
versus O3 0/4128 (0%).
