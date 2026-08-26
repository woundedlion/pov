# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_gsreactiondiffusion_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"`). Raw capture:
`build/prof/gsreactiondiffusion_o3.log`, captured 2026-08-26 01:25 local. This replaces `profile_gsreactiondiffusion_teensy_2026-08-09.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, flags `-D HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size: `FLASH: code:74,352, data:239,960, headers:8,248` / `RAM1: variables:315,232, code:40,808, padding:24,728, free:143,520` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 1185–1216
root counter cycles ÷ 600 MHz match the measured wall sum within **2.8 ppm**
(`tools/parse_profile.py build/prof/gsreactiondiffusion_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `grd_render` avg
48.85 ms/f, worst window 54.67 ms/f
(frames 1185–1216),
peak frame render 56.16 ms, spilled 0/2048
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 481–512)

```text
frame                          62.41ms  37.44Mcyc 100.0%
  grd_render                   54.48ms  32.69Mcyc  87.3%
    grd_rasterize              42.06ms  25.24Mcyc  67.4%
      grd_shader_draw          38.42ms  23.05Mcyc  61.6% x1.0 38422us/c
      grd_cull_flags            3.27ms   1.96Mcyc   5.2% x1.0 3271us/c
      grd_orient               372.5us  223.5kcyc   0.6% x1.0 372us/c
    grd_simulate               11.90ms   7.14Mcyc  19.1% x1.0 11905us/c
  rd_timeline_step              21.8us   13.1kcyc   0.0% x1.0 22us/c
  canvas_clear                  89.8us   53.9kcyc   0.1% x1.0 90us/c
  canvas_buffer_wait            7.82ms   4.69Mcyc  12.5% x1.0 7819us/c
```

Wall min/avg/max = 60.53/62.41/64.30 ms. The `grd_render` scope averages 54.48 ms/f in this window, while the exact frame-render peak is 56.16 ms. That is 6.34 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 1185–1216)

```text
frame                          62.48ms  37.49Mcyc 100.0%
  grd_render                   54.67ms  32.80Mcyc  87.5%
    grd_rasterize              42.77ms  25.66Mcyc  68.5%
      grd_shader_draw          38.99ms  23.39Mcyc  62.4% x1.0 38991us/c
      grd_cull_flags            3.40ms   2.04Mcyc   5.4% x1.0 3403us/c
      grd_orient               372.5us  223.5kcyc   0.6% x1.0 372us/c
    grd_simulate               11.90ms   7.14Mcyc  19.0% x1.0 11902us/c
  rd_timeline_step              24.2us   14.5kcyc   0.0% x1.0 24us/c
  canvas_clear                  89.6us   53.8kcyc   0.1% x1.0 90us/c
  canvas_buffer_wait            7.70ms   4.62Mcyc  12.3% x1.0 7696us/c
```

Wall min/avg/max = 60.57/62.48/64.49 ms. The `grd_render` scope averages 54.67 ms/f in this window, while the exact frame-render peak is 55.80 ms. That is 6.70 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.1/frame  min/avg/max 0.40/1.62/12.09 us  cpu 2.98%
isr_pack          143.9/frame  min/avg/max 6.08/7.32/10.35 us  cpu 1.68%
isr_dma_submit    143.9/frame  min/avg/max 0.69/0.94/3.75 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.87% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `grd_render` — 78.2% of measured root time, 48.85 ms/f average.
2. `grd_rasterize` — 59.0% of measured root time, 36.84 ms/f average.
3. `grd_shader_draw` — 53.0% of measured root time, 33.09 ms/f average.

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

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GSReactionDiffusion`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 55.26 ms versus
56.16 ms here: global O3 raises the peak by 0.90 ms (1.6%). O3-image minus shipping-image
size deltas are **FLASH code +13,224 B** and **ITCM code
+11,280 B**. Spill fractions are compared rather than raw counts:
shipping 0/2048 (0%)
versus O3 0/2048 (0%).
