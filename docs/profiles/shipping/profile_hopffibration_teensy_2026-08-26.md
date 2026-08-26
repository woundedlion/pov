# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HopfFibration`).
Raw capture: `build/prof/hopffibration_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_hopffibration_teensy_2026-07-30.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `HopfFibration::render_trails` and the Plot/raster `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh HopfFibration profile 70 32` |

Image size:

```text
FLASH: code:60,824, data:146,376, headers:8,864
       free for files:1,815,552
RAM1:  variables:315,200, code:40,904, padding:24,632
       free for local variables:143,552
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 193–224 root counter
cycles ÷ 600 MHz match the measured wall sum within **0.8 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hf_render_trails` averages
24.26 ms/frame; its worst window is 41.34 ms/frame
(frames 193–224). Peak frame render is
**48.50 ms** (frames 193–224);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 14.00 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 193–224)

```text
frame                         62.77 ms  37.66 Mcyc  100%
  hf_render_trails            41.34 ms  24.80 Mcyc   66%
    hf_trail_raster           32.26 ms  19.36 Mcyc   51%
      filter_blend             5.70 ms   3.42 Mcyc    9%  x49280.2 0.1us/c
    hf_trail_gate              8.24 ms   4.94 Mcyc   13%  x210.0 39.2us/c
  hf_project_record           186.3 us 111.80 kcyc    0%  x1.0 186.3us/c
  hf_advance_tumble             0.5 us     334 cyc    0%  x1.0 0.5us/c
  hf_timeline_step             17.3 us  10.41 kcyc    0%  x1.0 17.3us/c
  canvas_clear                 86.8 us  52.07 kcyc    0%  x1.0 86.8us/c
  canvas_buffer_wait          21.14 ms  12.69 Mcyc   34%  x1.0 21143.7us/c
```

Wall min/avg/max = 43.81/62.77/80.23 ms. `hf_render_trails`
accounts for 65.9% of this measured frame at 41.34 ms/frame.
The complete render is 41.63 ms; `canvas_buffer_wait` contributes
21.14 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–32)

```text
frame                         59.26 ms  35.56 Mcyc  100%
  hf_render_trails             9.68 ms   5.81 Mcyc   16%
    hf_trail_raster            6.82 ms   4.09 Mcyc   12%
      filter_blend             1.22 ms 731.44 kcyc    2%  x10523.9 0.1us/c
    hf_trail_gate              2.48 ms   1.49 Mcyc    4%  x203.4 12.2us/c
  hf_project_record           183.8 us 110.32 kcyc    0%  x1.0 183.8us/c
  hf_advance_tumble             0.5 us     325 cyc    0%  x1.0 0.5us/c
  hf_timeline_step             16.9 us  10.15 kcyc    0%  x1.0 16.9us/c
  canvas_clear                 85.9 us  51.54 kcyc    0%  x1.0 85.9us/c
  canvas_buffer_wait          49.29 ms  29.58 Mcyc   83%  x1.0 49293.2us/c
```

Wall min/avg/max = 0.31/59.26/69.16 ms. `hf_render_trails`
accounts for 16.3% of this measured frame at 9.68 ms/frame.
The complete render is 9.97 ms; `canvas_buffer_wait` contributes
49.29 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 49,280.2 pixels/frame, 4.75× the 10,368-pixel quadrant. `filter_blend` costs 69.3 cycles/blend; `hf_render_trails` uses 503.3 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.2/f  0.45/1.66/16.96 us  cpu 3.08%
isr_pack          143.8/f  6.24/6.93/9.81 us  cpu 1.59%
isr_dma_submit    143.8/f  0.59/0.94/2.27 us  cpu 0.22%
```

- Pack plus submit costs 1.13 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.89% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `hf_render_trails` — 38.9% of aggregated root time, 24.26 ms/frame: inclusive measured scope in the live driver.
2. `hf_trail_raster` — 27.6% of aggregated root time, 17.24 ms/frame: inclusive measured scope in the live driver.
3. `hf_trail_gate` — 9.9% of aggregated root time, 6.20 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches `HopfFibration::render_trails` and the Plot/raster `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HopfFibration` and
`HS_PROFILE_WINDOW=32`; `just profile HopfFibration` performs the
locked build, flash, capture, and artifact attestation.
