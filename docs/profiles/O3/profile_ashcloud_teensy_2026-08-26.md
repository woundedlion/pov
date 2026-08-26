# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_ashcloud_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh AshCloud profile_o3 70 32`). Raw capture:
`build/prof/ashcloud_o3.log`, captured 2026-08-26 02:45 local. This replaces `profile_ashcloud_teensy_2026-08-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AshCloud profile_o3 70 32` |

Image size: `FLASH: code:87,824, data:152,100, headers:8,908` / `RAM1: variables:315,104, code:38,248, padding:27,288, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 513–544
root counter cycles ÷ 600 MHz match the measured wall sum within **0.7 ppm**
(`tools/parse_profile.py build/prof/ashcloud_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
65.36 ms/f, worst window 70.91 ms/f
(frames 513–544),
peak frame render 79.81 ms, spilled 544/544
frames (100%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
The pass is red: 100% of frames exceed one display interval. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 513–544)

```text
frame                         124.94ms  74.97Mcyc 100.0%
  fx_shader_draw               70.91ms  42.55Mcyc  56.8% x1.0 70910us/c
  fx_prepare_frame              3.39ms   2.03Mcyc   2.7% x1.0 3388us/c
  fx_advance                    2.08ms   1.25Mcyc   1.7% x1.0 2081us/c
  fx_timeline_step              47.2us   28.3kcyc   0.0% x1.0 47us/c
  canvas_clear                  85.1us   51.0kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           48.42ms  29.05Mcyc  38.8% x1.0 48421us/c
```

Wall min/avg/max = 122.01/124.94/127.69 ms. The `fx_shader_draw` scope averages 70.91 ms/f in this window, while the exact frame-render peak is 79.81 ms. That is 17.31 ms of overrun against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         2304.7/frame  min/avg/max 0.39/5.12/23.33 us  cpu 9.43%
isr_pack          288.1/frame  min/avg/max 5.99/6.89/9.38 us  cpu 1.58%
isr_dma_submit    288.1/frame  min/avg/max 0.67/0.93/5.12 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 11.22% CPU. The worst render requires 27.7% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 52.4% of measured root time, 65.36 ms/f average.
2. `fx_prepare_frame` — 2.7% of measured root time, 3.37 ms/f average.
3. `fx_advance` — 1.7% of measured root time, 2.06 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AshCloud`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AshCloud profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 50.09 ms versus
79.81 ms here: global O3 raises the peak by 29.71 ms (59.3%). O3-image minus shipping-image
size deltas are **FLASH code +16,608 B** and **ITCM code
+12,112 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 544/544 (100%).
