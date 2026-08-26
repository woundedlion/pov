# KaleidoscopeStainedGlass on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopestainedglass_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeStainedGlass profile_o3 70 32`). Raw capture:
`build/prof/kaleidoscopestainedglass_o3.log`, captured 2026-08-26 02:51 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeStainedGlass 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeStainedGlass profile_o3 70 32` |

Image size: `FLASH: code:85,064, data:156,212, headers:8,580` / `RAM1: variables:315,104, code:37,672, padding:27,864, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 545–576
root counter cycles ÷ 600 MHz match the measured wall sum within **1.4 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopestainedglass_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
40.73 ms/f, worst window 43.02 ms/f
(frames 545–576),
peak frame render 46.99 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 545–576)

```text
frame                          62.44ms  37.46Mcyc 100.0%
  fx_shader_draw               43.02ms  25.81Mcyc  68.9% x1.0 43018us/c
  fx_prepare_frame              1.02ms  613.6kcyc   1.6% x1.0 1023us/c
  fx_advance                    2.25ms   1.35Mcyc   3.6% x1.0 2249us/c
  fx_timeline_step              35.9us   21.5kcyc   0.1% x1.0 36us/c
  canvas_clear                  87.7us   52.6kcyc   0.1% x1.0 88us/c
  canvas_buffer_wait           15.99ms   9.59Mcyc  25.6% x1.0 15991us/c
```

Wall min/avg/max = 61.51/62.44/63.28 ms. The `fx_shader_draw` scope averages 43.02 ms/f in this window, while the exact frame-render peak is 46.99 ms. That is 15.51 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.3/frame  min/avg/max 0.39/2.48/18.97 us  cpu 4.57%
isr_pack          144.0/frame  min/avg/max 6.39/7.22/9.59 us  cpu 1.66%
isr_dma_submit    144.0/frame  min/avg/max 0.62/0.93/7.02 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 6.44% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 65.3% of measured root time, 40.73 ms/f average.
2. `fx_advance` — 3.6% of measured root time, 2.23 ms/f average.
3. `fx_prepare_frame` — 1.5% of measured root time, 0.91 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `20ca3cb48892795a4575dd9d16d31a699f79df75`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeStainedGlass`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeStainedGlass profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 47.91 ms versus
46.99 ms here: global O3 lowers the peak by 0.93 ms (1.9%). O3-image minus shipping-image
size deltas are **FLASH code +13,784 B** and **ITCM code
+11,632 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
