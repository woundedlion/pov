# KaleidoscopeHexSoft on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopehexsoft_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopeHexSoft profile_o3 70 32`). Raw capture:
`build/prof/kaleidoscopehexsoft_o3.log`, captured 2026-08-26 02:26 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexSoft 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexSoft profile_o3 70 32` |

Image size: `FLASH: code:83,784, data:151,808, headers:9,144` / `RAM1: variables:315,104, code:37,848, padding:27,688, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 129–160
root counter cycles ÷ 600 MHz match the measured wall sum within **2.5 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopehexsoft_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
18.63 ms/f, worst window 18.70 ms/f
(frames 129–160),
peak frame render 24.58 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 961–992)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  fx_shader_draw               18.69ms  11.21Mcyc  29.9% x1.0 18687us/c
  fx_prepare_frame              3.32ms   1.99Mcyc   5.3% x1.0 3319us/c
  fx_advance                    2.24ms   1.35Mcyc   3.6% x1.0 2244us/c
  fx_timeline_step              37.3us   22.4kcyc   0.1% x1.0 37us/c
  canvas_clear                  85.9us   51.6kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           38.05ms  22.83Mcyc  60.9% x1.0 38047us/c
```

Wall min/avg/max = 62.09/62.43/62.69 ms. The `fx_shader_draw` scope averages 18.69 ms/f in this window, while the exact frame-render peak is 24.58 ms. That is 37.92 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.0/frame  min/avg/max 0.39/1.51/22.94 us  cpu 2.77%
isr_pack          144.0/frame  min/avg/max 5.99/6.83/9.15 us  cpu 1.57%
isr_dma_submit    144.0/frame  min/avg/max 0.65/0.94/8.46 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.55% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 29.9% of measured root time, 18.63 ms/f average.
2. `fx_prepare_frame` — 5.4% of measured root time, 3.35 ms/f average.
3. `fx_advance` — 3.5% of measured root time, 2.21 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexSoft`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopeHexSoft profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 27.30 ms versus
24.58 ms here: global O3 lowers the peak by 2.72 ms (10.0%). O3-image minus shipping-image
size deltas are **FLASH code +14,496 B** and **ITCM code
+11,712 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
