# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_gnomonicstars_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh GnomonicStars profile_o3 70 32`). Raw capture:
`build/prof/gnomonicstars_o3.log`, captured 2026-08-26 01:22 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GnomonicStars profile_o3 70 32` |

Image size: `FLASH: code:63,160, data:148,732, headers:8,268` / `RAM1: variables:315,136, code:40,568, padding:24,968, free:143,616` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 193–224
root counter cycles ÷ 600 MHz match the measured wall sum within **1.2 ppm**
(`tools/parse_profile.py build/prof/gnomonicstars_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `gn_draw_stars` avg
8.29 ms/f, worst window 19.05 ms/f
(frames 193–224),
peak frame render 26.29 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 193–224)

```text
frame                          62.40ms  37.44Mcyc 100.0%
  gn_draw_stars                19.05ms  11.43Mcyc  30.5%
    gn_star_scan               18.55ms  11.13Mcyc  29.7%
      filter_blend             811.1us  486.7kcyc   1.3% x9118.5 0us/c
  gn_timeline_step              37.7us   22.6kcyc   0.1% x1.0 38us/c
  canvas_clear                  85.0us   51.0kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           43.23ms  25.94Mcyc  69.3% x1.0 43227us/c
```

Wall min/avg/max = 52.00/62.40/72.14 ms. The `gn_draw_stars` scope averages 19.05 ms/f in this window, while the exact frame-render peak is 26.29 ms. That is 36.21 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 9,118.5 blended px/frame versus the 10,368-px quadrant (0.88× coverage), at 53.4 cyc/blend. `gn_draw_stars` contributes 1253.6 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.2/frame  min/avg/max 0.47/1.55/13.35 us  cpu 2.86%
isr_pack          143.9/frame  min/avg/max 6.03/6.68/9.09 us  cpu 1.53%
isr_dma_submit    143.9/frame  min/avg/max 0.77/0.93/6.03 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.60% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `gn_draw_stars` — 13.3% of measured root time, 8.29 ms/f average.
2. `gn_star_scan` — 12.5% of measured root time, 7.79 ms/f average.
3. `filter_blend` — 0.5% of measured root time, 0.29 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `gn_star_scan`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GnomonicStars`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GnomonicStars profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 29.64 ms versus
26.29 ms here: global O3 lowers the peak by 3.35 ms (11.3%). O3-image minus shipping-image
size deltas are **FLASH code +12,184 B** and **ITCM code
+11,104 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
