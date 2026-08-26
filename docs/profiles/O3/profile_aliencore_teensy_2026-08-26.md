# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_aliencore_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh AlienCore profile_o3 70 32`). Raw capture:
`build/prof/aliencore_o3.log`, captured 2026-08-26 02:30 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienCore 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienCore profile_o3 70 32` |

Image size: `FLASH: code:82,952, data:151,772, headers:8,988` / `RAM1: variables:315,104, code:37,608, padding:27,928, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 865–896
root counter cycles ÷ 600 MHz match the measured wall sum within **2.7 ppm**
(`tools/parse_profile.py build/prof/aliencore_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
15.15 ms/f, worst window 15.16 ms/f
(frames 865–896),
peak frame render 20.10 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 1–32)

```text
frame                          59.71ms  35.83Mcyc 100.0%
  fx_shader_draw               15.15ms   9.09Mcyc  25.4% x1.0 15147us/c
  fx_prepare_frame              89.8us   53.9kcyc   0.2% x1.0 90us/c
  fx_advance                    2.35ms   1.41Mcyc   3.9% x1.0 2354us/c
  fx_timeline_step              34.9us   20.9kcyc   0.1% x1.0 35us/c
  canvas_clear                  86.2us   51.7kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           41.99ms  25.19Mcyc  70.3% x1.0 41988us/c
```

Wall min/avg/max = 17.57/59.71/62.60 ms. The `fx_shader_draw` scope averages 15.15 ms/f in this window, while the exact frame-render peak is 20.10 ms. That is 42.40 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1378.4/frame  min/avg/max 0.40/1.34/18.82 us  cpu 2.48%
isr_pack          136.3/frame  min/avg/max 0.52/7.10/9.50 us  cpu 1.29%
isr_dma_submit    136.3/frame  min/avg/max 0.63/0.93/7.81 us  cpu 0.17%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 3.94% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 24.3% of measured root time, 15.15 ms/f average.
2. `fx_advance` — 3.6% of measured root time, 2.25 ms/f average.
3. `canvas_clear` — 0.1% of measured root time, 0.09 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienCore`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AlienCore profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 21.09 ms versus
20.10 ms here: global O3 lowers the peak by 0.99 ms (4.7%). O3-image minus shipping-image
size deltas are **FLASH code +13,648 B** and **ITCM code
+11,440 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
