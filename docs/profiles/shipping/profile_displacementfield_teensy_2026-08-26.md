# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/prof/displacementfield_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_displacementfield_teensy_2026-08-18.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `DisplacementField::draw_rings`, field/hue helpers, and ring-SDF `HS_O3_FN` paths |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile 70 32` |

Image size:

```text
FLASH: code:78,496, data:150,388, headers:8,684
       free for files:1,794,048
RAM1:  variables:315,232, code:43,464, padding:22,072
       free for local variables:143,520
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 609–640 root counter
cycles ÷ 600 MHz match the measured wall sum within **3.2 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `df_timeline_step` averages
46.83 ms/frame; its worst window is 54.99 ms/frame
(frames 609–640). Peak frame render is
**57.40 ms** (frames 577–608);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 5.10 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 609–640)

```text
frame                         62.46 ms  37.48 Mcyc  100%
  df_timeline_step            54.99 ms  32.99 Mcyc   88%
    df_draw_rings             54.89 ms  32.93 Mcyc   88%
      df_hue_table_prep        1.78 ms   1.07 Mcyc    3%  x21.5 82.6us/c
      df_lut_bake             10.30 ms   6.18 Mcyc   16%  x43.2 238.3us/c
      df_chunk_cull            1.05 ms 632.84 kcyc    2%  x45.4 23.2us/c
      df_fused_scan           40.30 ms  24.18 Mcyc   65%
        filter_blend           1.17 ms 701.54 kcyc    2%  x10109.4 0.1us/c
  canvas_clear                 89.7 us  53.81 kcyc    0%  x1.0 89.7us/c
  canvas_buffer_wait           7.38 ms   4.43 Mcyc   12%  x1.0 7383.3us/c
  df_prepare_fields             0.2 us     163 cyc    0%  x1.0 0.2us/c
```

Wall min/avg/max = 60.46/62.46/64.10 ms. `df_timeline_step`
accounts for 88.0% of this measured frame at 54.99 ms/frame.
The complete render is 55.08 ms; `canvas_buffer_wait` contributes
7.38 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 897–928)

```text
frame                         62.94 ms  37.76 Mcyc  100%
  df_timeline_step            20.55 ms  12.33 Mcyc   33%
    df_draw_rings             20.43 ms  12.26 Mcyc   32%
      df_hue_table_prep       464.9 us 278.95 kcyc    1%  x11.2 41.3us/c
      df_lut_bake              2.07 ms   1.24 Mcyc    3%  x13.7 151.5us/c
      df_chunk_cull           371.5 us 222.90 kcyc    1%  x14.5 25.7us/c
      df_fused_scan           16.43 ms   9.86 Mcyc   26%
        filter_blend           1.02 ms 611.70 kcyc    2%  x9145.9 0.1us/c
  canvas_clear                 85.4 us  51.26 kcyc    0%  x1.0 85.4us/c
  canvas_buffer_wait          42.29 ms  25.37 Mcyc   67%  x1.0 42290.0us/c
  df_prepare_fields             1.3 us     807 cyc    0%  x1.0 1.3us/c
```

Wall min/avg/max = 57.91/62.94/64.88 ms. `df_timeline_step`
accounts for 32.7% of this measured frame at 20.55 ms/frame.
The complete render is 20.65 ms; `canvas_buffer_wait` contributes
42.29 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 10,109.4 pixels/frame, 0.98× the 10,368-pixel quadrant. `filter_blend` costs 69.4 cycles/blend; `df_timeline_step` uses 3263.5 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.1/f  0.43/1.68/20.42 us  cpu 3.11%
isr_pack          143.8/f  6.24/6.88/17.63 us  cpu 1.58%
isr_dma_submit    143.8/f  0.60/0.94/8.00 us  cpu 0.22%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.91% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `df_timeline_step` — 75.1% of aggregated root time, 46.83 ms/frame: inclusive measured scope in the live driver.
2. `df_draw_rings` — 74.9% of aggregated root time, 46.72 ms/frame: inclusive measured scope in the live driver.
3. `df_fused_scan` — 54.0% of aggregated root time, 33.66 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches `DisplacementField::draw_rings`, field/hue helpers, and ring-SDF `HS_O3_FN` paths; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=DisplacementField` and
`HS_PROFILE_WINDOW=32`; `just profile DisplacementField` performs the
locked build, flash, capture, and artifact attestation.
