# BZReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_bzreactiondiffusion_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh BZReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"`). Raw capture:
`build/prof/bzreactiondiffusion_o3.log`, captured 2026-08-26 01:16 local. This replaces `profile_bzreactiondiffusion_teensy_2026-08-03.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | BZReactionDiffusion 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, flags `-D HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh BZReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:74,224, data:239,776, headers:8,560` / `RAM1: variables:315,136, code:39,640, padding:25,896, free:143,616` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 1089–1120
root counter cycles ÷ 600 MHz match the measured wall sum within **1.0 ppm**
(`tools/parse_profile.py build/prof/bzreactiondiffusion_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `bz_render` avg
46.93 ms/f, worst window 47.69 ms/f
(frames 1089–1120),
peak frame render 48.65 ms, spilled 0/2048
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 1441–1472)

```text
frame                          62.48ms  37.49Mcyc 100.0%
  bz_render                    47.66ms  28.59Mcyc  76.3%
    bz_raster                  42.83ms  25.70Mcyc  68.5% x1.0 42826us/c
    bz_orient                  372.4us  223.4kcyc   0.6% x1.0 372us/c
    bz_physics                  4.46ms   2.67Mcyc   7.1% x1.0 4458us/c
  rd_timeline_step              34.3us   20.6kcyc   0.1% x1.0 34us/c
  canvas_clear                  87.8us   52.7kcyc   0.1% x1.0 88us/c
  canvas_buffer_wait           14.70ms   8.82Mcyc  23.5% x1.0 14697us/c
```

Wall min/avg/max = 61.07/62.48/63.95 ms. The `bz_render` scope averages 47.66 ms/f in this window, while the exact frame-render peak is 48.65 ms. That is 13.85 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.5/frame  min/avg/max 0.33/1.57/13.54 us  cpu 2.89%
isr_pack          144.1/frame  min/avg/max 6.13/6.96/9.94 us  cpu 1.60%
isr_dma_submit    144.1/frame  min/avg/max 0.66/0.93/5.94 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.70% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `bz_render` — 75.2% of measured root time, 46.93 ms/f average.
2. `bz_raster` — 67.4% of measured root time, 42.10 ms/f average.
3. `bz_physics` — 7.1% of measured root time, 4.46 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image activates effect-local `HS_O3_FN` simulation/raster regions plus any shared hot paths; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=1200` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=BZReactionDiffusion`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh BZReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 48.91 ms versus
48.65 ms here: global O3 lowers the peak by 0.27 ms (0.5%). O3-image minus shipping-image
size deltas are **FLASH code +12,760 B** and **ITCM code
+10,224 B**. Spill fractions are compared rather than raw counts:
shipping 0/2048 (0%)
versus O3 0/2048 (0%).
