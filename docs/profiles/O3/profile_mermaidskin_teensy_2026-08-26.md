# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_mermaidskin_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh MermaidSkin profile_o3 70 32`). Raw capture:
`build/prof/mermaidskin_o3.log`, captured 2026-08-26 02:43 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MermaidSkin 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh MermaidSkin profile_o3 70 32` |

Image size: `FLASH: code:87,960, data:151,932, headers:8,940` / `RAM1: variables:315,104, code:38,088, padding:27,448, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 641–672
root counter cycles ÷ 600 MHz match the measured wall sum within **0.8 ppm**
(`tools/parse_profile.py build/prof/mermaidskin_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
47.32 ms/f, worst window 47.57 ms/f
(frames 641–672),
peak frame render 54.55 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 289–320)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  fx_shader_draw               47.56ms  28.54Mcyc  76.2% x1.0 47565us/c
  fx_prepare_frame              4.05ms   2.43Mcyc   6.5% x1.0 4053us/c
  fx_advance                    2.00ms   1.20Mcyc   3.2% x1.0 2005us/c
  fx_timeline_step              36.6us   22.0kcyc   0.1% x1.0 37us/c
  canvas_clear                  88.9us   53.3kcyc   0.1% x1.0 89us/c
  canvas_buffer_wait            8.65ms   5.19Mcyc  13.9% x1.0 8652us/c
```

Wall min/avg/max = 61.29/62.42/64.05 ms. The `fx_shader_draw` scope averages 47.56 ms/f in this window, while the exact frame-render peak is 54.55 ms. That is 7.95 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame  min/avg/max 0.39/4.65/21.53 us  cpu 8.56%
isr_pack          144.0/frame  min/avg/max 6.47/7.42/9.49 us  cpu 1.71%
isr_dma_submit    144.0/frame  min/avg/max 0.69/0.93/4.50 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 10.48% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 75.8% of measured root time, 47.32 ms/f average.
2. `fx_prepare_frame` — 5.9% of measured root time, 3.66 ms/f average.
3. `fx_advance` — 3.2% of measured root time, 2.01 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MermaidSkin`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh MermaidSkin profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 45.60 ms versus
54.55 ms here: global O3 raises the peak by 8.95 ms (19.6%). O3-image minus shipping-image
size deltas are **FLASH code +16,544 B** and **ITCM code
+11,952 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
