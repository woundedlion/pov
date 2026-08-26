# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_cosmiceyeball_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh CosmicEyeball profile_o3 70 32`). Raw capture:
`build/prof/cosmiceyeball_o3.log`, captured 2026-08-26 03:05 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | CosmicEyeball 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh CosmicEyeball profile_o3 70 32` |

Image size: `FLASH: code:80,984, data:151,744, headers:8,936` / `RAM1: variables:315,104, code:36,296, padding:29,240, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 545–576
root counter cycles ÷ 600 MHz match the measured wall sum within **3.0 ppm**
(`tools/parse_profile.py build/prof/cosmiceyeball_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
20.02 ms/f, worst window 20.04 ms/f
(frames 545–576),
peak frame render 22.97 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 577–608)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  fx_shader_draw               20.03ms  12.02Mcyc  32.1% x1.0 20028us/c
  fx_prepare_frame             370.9us  222.5kcyc   0.6% x1.0 371us/c
  fx_advance                    2.20ms   1.32Mcyc   3.5% x1.0 2198us/c
  fx_timeline_step              29.1us   17.5kcyc   0.0% x1.0 29us/c
  canvas_clear                  85.7us   51.4kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           39.70ms  23.82Mcyc  63.6% x1.0 39703us/c
```

Wall min/avg/max = 59.78/62.42/65.03 ms. The `fx_shader_draw` scope averages 20.03 ms/f in this window, while the exact frame-render peak is 22.97 ms. That is 39.53 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame  min/avg/max 0.39/1.53/26.23 us  cpu 2.82%
isr_pack          144.0/frame  min/avg/max 6.17/7.04/9.77 us  cpu 1.62%
isr_dma_submit    144.0/frame  min/avg/max 0.67/0.94/9.13 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.65% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 32.1% of measured root time, 20.02 ms/f average.
2. `fx_advance` — 3.6% of measured root time, 2.23 ms/f average.
3. `fx_prepare_frame` — 0.6% of measured root time, 0.38 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=CosmicEyeball`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh CosmicEyeball profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 24.01 ms versus
22.97 ms here: global O3 lowers the peak by 1.05 ms (4.4%). O3-image minus shipping-image
size deltas are **FLASH code +12,552 B** and **ITCM code
+10,320 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
