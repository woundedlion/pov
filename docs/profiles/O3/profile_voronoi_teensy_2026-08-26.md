# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_voronoi_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh Voronoi profile_o3 70 32`). Raw capture:
`build/prof/voronoi_o3.log`, captured 2026-08-26 01:39 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Voronoi 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Voronoi profile_o3 70 32` |

Image size: `FLASH: code:54,912, data:145,952, headers:9,056` / `RAM1: variables:315,040, code:33,160, padding:32,376, free:143,712` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 353–384
root counter cycles ÷ 600 MHz match the measured wall sum within **0.2 ppm**
(`tools/parse_profile.py build/prof/voronoi_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `vo_shade` avg
7.10 ms/f, worst window 7.18 ms/f
(frames 353–384),
peak frame render 7.71 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 769–800)

```text
frame                          62.45ms  37.47Mcyc 100.0%
  vo_shade                      7.17ms   4.30Mcyc  11.5% x1.0 7173us/c
  vo_kdtree                    258.2us  154.9kcyc   0.4% x1.0 258us/c
  vo_animate                    40.7us   24.4kcyc   0.1% x1.0 41us/c
  canvas_clear                  84.9us   51.0kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           54.89ms  32.93Mcyc  87.9% x1.0 54888us/c
```

Wall min/avg/max = 62.26/62.45/62.65 ms. The `vo_shade` scope averages 7.17 ms/f in this window, while the exact frame-render peak is 7.71 ms. That is 54.79 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame  min/avg/max 0.40/1.50/11.96 us  cpu 2.76%
isr_pack          144.0/frame  min/avg/max 6.43/7.00/9.44 us  cpu 1.61%
isr_dma_submit    144.0/frame  min/avg/max 0.62/0.94/1.13 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.58% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `vo_shade` — 11.4% of measured root time, 7.10 ms/f average.
2. `vo_kdtree` — 0.4% of measured root time, 0.26 ms/f average.
3. `canvas_clear` — 0.1% of measured root time, 0.08 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- No effect-local selective-O3 boundary is exposed by this counter tree; the paired selective-O3 capture is the shipping reference, while global O3 compiles every translation unit.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Voronoi`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Voronoi profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 8.96 ms versus
7.71 ms here: global O3 lowers the peak by 1.25 ms (13.9%). O3-image minus shipping-image
size deltas are **FLASH code +15,568 B** and **ITCM code
+12,688 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
