# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_displacementfield_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh DisplacementField profile_o3 70 32`). Raw capture:
`build/prof/displacementfield_o3.log`, captured 2026-08-26 01:20 local. This replaces `profile_displacementfield_teensy_2026-08-18.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile_o3 70 32` |

Image size: `FLASH: code:104,008, data:150,096, headers:9,064` / `RAM1: variables:315,264, code:65,384, padding:152, free:143,488` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 737–768
root counter cycles ÷ 600 MHz match the measured wall sum within **0.1 ppm**
(`tools/parse_profile.py build/prof/displacementfield_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `df_timeline_step` avg
46.71 ms/f, worst window 54.69 ms/f
(frames 737–768),
peak frame render 57.48 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 577–608)

```text
frame                          62.55ms  37.53Mcyc 100.0%
  df_timeline_step             53.45ms  32.07Mcyc  85.5%
    df_draw_rings              53.37ms  32.02Mcyc  85.3%
      df_hue_table_prep         1.64ms  983.0kcyc   2.6% x22.1 74us/c
      df_lut_bake               9.64ms   5.78Mcyc  15.4% x35.7 270us/c
      df_chunk_cull            950.7us  570.4kcyc   1.5% x41.0 23us/c
      df_fused_scan            39.88ms  23.93Mcyc  63.8%
        filter_blend            1.06ms  635.7kcyc   1.7% x9951.0 0us/c
  canvas_clear                  88.6us   53.1kcyc   0.1% x1.0 89us/c
  canvas_buffer_wait            9.00ms   5.40Mcyc  14.4% x1.0 9003us/c
  df_prepare_fields              0.2us     163cyc   0.0% x1.0 0us/c
```

Wall min/avg/max = 55.73/62.55/69.24 ms. The `df_timeline_step` scope averages 53.45 ms/f in this window, while the exact frame-render peak is 57.48 ms. That is 5.02 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 737–768)

```text
frame                          62.45ms  37.47Mcyc 100.0%
  df_timeline_step             54.69ms  32.81Mcyc  87.6%
    df_draw_rings              54.60ms  32.76Mcyc  87.4%
      df_hue_table_prep         1.46ms  876.5kcyc   2.3% x21.3 69us/c
      df_lut_bake              10.33ms   6.20Mcyc  16.5% x46.7 221us/c
      df_chunk_cull             1.07ms  639.7kcyc   1.7% x47.3 23us/c
      df_fused_scan            40.30ms  24.18Mcyc  64.5%
        filter_blend            1.08ms  649.2kcyc   1.7% x10165.2 0us/c
  canvas_clear                  89.4us   53.7kcyc   0.1% x1.0 89us/c
  canvas_buffer_wait            7.68ms   4.61Mcyc  12.3% x1.0 7677us/c
  df_prepare_fields              0.3us     180cyc   0.0% x1.0 0us/c
```

Wall min/avg/max = 60.41/62.45/63.70 ms. The `df_timeline_step` scope averages 54.69 ms/f in this window, while the exact frame-render peak is 56.20 ms. That is 6.30 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 9,951.0 blended px/frame versus the 10,368-px quadrant (0.96× coverage), at 63.9 cyc/blend. `df_timeline_step` contributes 3223.0 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1153.7/frame  min/avg/max 0.35/1.59/11.82 us  cpu 2.93%
isr_pack          144.2/frame  min/avg/max 6.22/6.95/9.31 us  cpu 1.60%
isr_dma_submit    144.2/frame  min/avg/max 0.65/0.94/1.02 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.74% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `df_timeline_step` — 74.9% of measured root time, 46.71 ms/f average.
2. `df_draw_rings` — 74.7% of measured root time, 46.61 ms/f average.
3. `df_fused_scan` — 54.6% of measured root time, 34.07 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `df_fused_scan`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image activates effect-local `HS_O3_FN` simulation/raster regions plus any shared hot paths; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh DisplacementField profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 57.40 ms versus
57.48 ms here: global O3 raises the peak by 0.08 ms (0.1%). O3-image minus shipping-image
size deltas are **FLASH code +25,512 B** and **ITCM code
+21,920 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
