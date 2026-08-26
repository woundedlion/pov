# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopesmooth_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeSmooth profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"`). Raw capture:
`build/prof/kaleidoscopesmooth_o3.log`, captured 2026-08-26 02:56 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, flags `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeSmooth profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size: `FLASH: code:84,528, data:151,984, headers:8,224` / `RAM1: variables:315,104, code:37,848, padding:27,688, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 881–896
root counter cycles ÷ 600 MHz match the measured wall sum within **0.9 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopesmooth_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
24.26 ms/f, worst window 29.07 ms/f
(frames 881–896),
peak frame render 32.58 ms, spilled 0/4128
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 881–896)

```text
frame                          62.41ms  37.45Mcyc 100.0%
  fx_shader_draw               29.07ms  17.44Mcyc  46.6% x1.0 29073us/c
  fx_prepare_frame             714.2us  428.5kcyc   1.1% x1.0 714us/c
  fx_advance                    2.06ms   1.24Mcyc   3.3% x1.0 2061us/c
  fx_timeline_step              94.0us   56.4kcyc   0.2% x1.0 94us/c
  canvas_clear                  86.1us   51.6kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           30.37ms  18.22Mcyc  48.7% x1.0 30374us/c
```

Wall min/avg/max = 61.31/62.41/63.34 ms. The `fx_shader_draw` scope averages 29.07 ms/f in this window, while the exact frame-render peak is 32.58 ms. That is 29.92 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 3169–3184)

```text
frame                          62.47ms  37.48Mcyc 100.0%
  fx_shader_draw               21.26ms  12.76Mcyc  34.0% x1.0 21264us/c
  fx_prepare_frame             908.4us  545.0kcyc   1.5% x1.0 908us/c
  fx_advance                    2.05ms   1.23Mcyc   3.3% x1.0 2048us/c
  fx_timeline_step             102.8us   61.7kcyc   0.2% x1.0 103us/c
  canvas_clear                  85.6us   51.3kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           38.04ms  22.82Mcyc  60.9% x1.0 38041us/c
```

Wall min/avg/max = 61.66/62.47/63.34 ms. The `fx_shader_draw` scope averages 21.26 ms/f in this window, while the exact frame-render peak is 24.92 ms. That is 37.58 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `fx_shader_draw` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `2` | — | 67/67 | — | 29.07 | 32.58 | 16.0 |
| `1` | — | 18/18 | — | 28.60 | 32.26 | 16.0 |
| `0` | — | 38/38 | — | 27.66 | 30.49 | 16.0 |
| `4` | — | 67/67 | — | 25.16 | 28.56 | 16.0 |
| `3` | — | 68/68 | — | 25.10 | 28.89 | 16.0 |

220 advance markers visit 4/4 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.9/frame  min/avg/max 0.39/1.50/10.64 us  cpu 2.76%
isr_pack          143.9/frame  min/avg/max 5.99/6.64/8.80 us  cpu 1.52%
isr_dma_submit    143.9/frame  min/avg/max 0.59/0.93/3.53 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.49% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 38.9% of measured root time, 24.26 ms/f average.
2. `fx_advance` — 3.2% of measured root time, 2.01 ms/f average.
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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeSmooth`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh KaleidoscopeSmooth profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 35.85 ms versus
32.58 ms here: global O3 lowers the peak by 3.27 ms (9.1%). O3-image minus shipping-image
size deltas are **FLASH code +14,560 B** and **ITCM code
+11,712 B**. Spill fractions are compared rather than raw counts:
shipping 0/4128 (0%)
versus O3 0/4128 (0%).
