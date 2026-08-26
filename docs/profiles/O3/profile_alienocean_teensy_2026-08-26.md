# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_alienocean_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh AlienOcean profile_o3 70 32`). Raw capture:
`build/prof/alienocean_o3.log`, captured 2026-08-26 02:28 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienOcean 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienOcean profile_o3 70 32` |

Image size: `FLASH: code:83,032, data:151,772, headers:8,908` / `RAM1: variables:315,104, code:37,608, padding:27,928, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 129–160
root counter cycles ÷ 600 MHz match the measured wall sum within **1.7 ppm**
(`tools/parse_profile.py build/prof/alienocean_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
19.83 ms/f, worst window 19.91 ms/f
(frames 129–160),
peak frame render 24.98 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 1–32)

```text
frame                          60.02ms  36.01Mcyc 100.0%
  fx_shader_draw               19.82ms  11.89Mcyc  33.0% x1.0 19815us/c
  fx_prepare_frame             467.4us  280.4kcyc   0.8% x1.0 467us/c
  fx_advance                    2.22ms   1.33Mcyc   3.7% x1.0 2221us/c
  fx_timeline_step              34.7us   20.8kcyc   0.1% x1.0 35us/c
  canvas_clear                  86.4us   51.9kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           37.38ms  22.43Mcyc  62.3% x1.0 37380us/c
```

Wall min/avg/max = 22.43/60.02/62.64 ms. The `fx_shader_draw` scope averages 19.82 ms/f in this window, while the exact frame-render peak is 24.98 ms. That is 37.52 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1381.2/frame  min/avg/max 0.40/1.34/17.90 us  cpu 2.47%
isr_pack          136.7/frame  min/avg/max 0.52/6.97/9.44 us  cpu 1.27%
isr_dma_submit    136.7/frame  min/avg/max 0.62/0.93/3.81 us  cpu 0.16%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 3.90% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 31.8% of measured root time, 19.83 ms/f average.
2. `fx_advance` — 3.5% of measured root time, 2.20 ms/f average.
3. `fx_prepare_frame` — 0.7% of measured root time, 0.44 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AlienOcean`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh AlienOcean profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 27.79 ms versus
24.98 ms here: global O3 lowers the peak by 2.81 ms (10.1%). O3-image minus shipping-image
size deltas are **FLASH code +13,648 B** and **ITCM code
+11,440 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
