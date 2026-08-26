# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_alienbrain_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh AlienBrain profile_o3 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"`). Raw capture:
`build/prof/alienbrain_o3.log`, captured 2026-08-26 02:24 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienBrain 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 300 s capture, flags `-D HS_PROFILE_EPOCH_REVS=2600` |
| Reproduce | `bash tools/profile_one.sh AlienBrain profile_o3 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` |

Image size: `FLASH: code:84,336, data:151,904, headers:8,496` / `RAM1: variables:315,104, code:37,848, padding:27,688, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 4561–4576
root counter cycles ÷ 600 MHz match the measured wall sum within **3.3 ppm**
(`tools/parse_profile.py build/prof/alienbrain_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
22.18 ms/f, worst window 22.83 ms/f
(frames 4561–4576),
peak frame render 27.90 ms, spilled 0/4768
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 1–16)

```text
frame                          57.96ms  34.78Mcyc 100.0%
  fx_shader_draw               22.55ms  13.53Mcyc  38.9% x1.0 22549us/c
  fx_prepare_frame             767.3us  460.4kcyc   1.3% x1.0 767us/c
  fx_advance                    2.27ms   1.36Mcyc   3.9% x1.0 2273us/c
  fx_timeline_step              38.8us   23.3kcyc   0.1% x1.0 39us/c
  canvas_clear                  86.9us   52.2kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           32.24ms  19.34Mcyc  55.6% x1.0 32236us/c
```

Wall min/avg/max = 25.55/57.96/62.55 ms. The `fx_shader_draw` scope averages 22.55 ms/f in this window, while the exact frame-render peak is 27.90 ms. That is 34.60 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 2721–2736)

```text
frame                          62.42ms  37.45Mcyc 100.0%
  fx_shader_draw               20.80ms  12.48Mcyc  33.3% x1.0 20801us/c
  fx_prepare_frame             571.6us  343.0kcyc   0.9% x1.0 572us/c
  fx_advance                    2.23ms   1.34Mcyc   3.6% x1.0 2234us/c
  fx_timeline_step              41.4us   24.9kcyc   0.1% x1.0 41us/c
  canvas_clear                  85.4us   51.3kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           38.68ms  23.21Mcyc  62.0% x1.0 38677us/c
```

Wall min/avg/max = 62.20/62.42/62.56 ms. The `fx_shader_draw` scope averages 20.80 ms/f in this window, while the exact frame-render peak is 23.82 ms. That is 38.67 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `1` | — | 58/58 | — | 22.83 | 26.09 | 16.0 |
| `2` | — | 67/67 | — | 22.81 | 26.07 | 16.0 |
| `0` | — | 38/38 | — | 22.81 | 27.90 | 16.0 |
| `4` | — | 67/67 | — | 22.76 | 26.16 | 16.0 |
| `3` | — | 68/68 | — | 22.04 | 25.73 | 16.0 |

260 advance markers visit 4/4 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1613.7/frame  min/avg/max 0.33/1.11/18.42 us  cpu 2.05%
isr_pack          129.8/frame  min/avg/max 0.52/6.57/8.88 us  cpu 0.97%
isr_dma_submit    129.8/frame  min/avg/max 0.62/0.94/9.08 us  cpu 0.13%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 3.15% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 35.6% of measured root time, 22.18 ms/f average.
2. `fx_advance` — 3.6% of measured root time, 2.24 ms/f average.
3. `fx_prepare_frame` — 1.1% of measured root time, 0.69 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=2600` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienBrain`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh AlienBrain profile_o3 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 31.61 ms versus
27.90 ms here: global O3 lowers the peak by 3.71 ms (11.7%). O3-image minus shipping-image
size deltas are **FLASH code +14,560 B** and **ITCM code
+11,712 B**. Spill fractions are compared rather than raw counts:
shipping 0/4768 (0%)
versus O3 0/4768 (0%).
