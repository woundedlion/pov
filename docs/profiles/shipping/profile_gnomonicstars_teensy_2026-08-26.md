# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GnomonicStars`).
Raw capture: `build/prof/gnomonicstars_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_gnomonicstars_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses shape scan, transformer, and raster `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GnomonicStars profile 70 32` |

Image size:

```text
FLASH: code:50,976, data:148,756, headers:9,164
       free for files:1,822,720
RAM1:  variables:315,136, code:29,464, padding:3,304
       free for local variables:176,384
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 193–224 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.3 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `gn_draw_stars` averages
10.12 ms/frame; its worst window is 21.72 ms/frame
(frames 193–224). Peak frame render is
**29.64 ms** (frames 193–224);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 32.86 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 193–224)

```text
frame                         62.41 ms  37.44 Mcyc  100%
  gn_draw_stars               21.72 ms  13.03 Mcyc   35%
    gn_star_scan              21.05 ms  12.63 Mcyc   34%
      filter_blend            909.3 us 545.59 kcyc    1%  x9118.5 0.1us/c
  gn_timeline_step             38.2 us  22.92 kcyc    0%  x1.0 38.2us/c
  canvas_clear                 85.3 us  51.19 kcyc    0%  x1.0 85.3us/c
  canvas_buffer_wait          40.56 ms  24.34 Mcyc   65%  x1.0 40562.4us/c
```

Wall min/avg/max = 50.89/62.41/73.19 ms. `gn_draw_stars`
accounts for 34.8% of this measured frame at 21.72 ms/frame.
The complete render is 21.84 ms; `canvas_buffer_wait` contributes
40.56 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 417–448)

```text
frame                         62.46 ms  37.47 Mcyc  100%
  gn_draw_stars                8.58 ms   5.15 Mcyc   14%
    gn_star_scan               7.91 ms   4.75 Mcyc   13%
      filter_blend            244.9 us 146.96 kcyc    0%  x2308.9 0.1us/c
  gn_timeline_step             37.1 us  22.25 kcyc    0%  x1.0 37.1us/c
  canvas_clear                 84.9 us  50.95 kcyc    0%  x1.0 84.9us/c
  canvas_buffer_wait          53.75 ms  32.25 Mcyc   86%  x1.0 53752.3us/c
```

Wall min/avg/max = 58.45/62.46/66.50 ms. `gn_draw_stars`
accounts for 13.7% of this measured frame at 8.58 ms/frame.
The complete render is 8.70 ms; `canvas_buffer_wait` contributes
53.75 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 9,118.5 pixels/frame, 0.88× the 10,368-pixel quadrant. `filter_blend` costs 59.8 cycles/blend; `gn_draw_stars` uses 1429.2 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.6/f  0.50/1.64/20.07 us  cpu 3.03%
isr_pack          143.8/f  6.24/6.79/9.64 us  cpu 1.56%
isr_dma_submit    143.8/f  0.62/0.94/2.08 us  cpu 0.22%
```

- Pack plus submit costs 1.11 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.81% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `gn_draw_stars` — 16.2% of aggregated root time, 10.12 ms/frame: inclusive measured scope in the live driver.
2. `gn_star_scan` — 15.1% of aggregated root time, 9.45 ms/frame: inclusive measured scope in the live driver.
3. `filter_blend` — 0.5% of aggregated root time, 0.33 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches shape scan, transformer, and raster `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=GnomonicStars` and
`HS_PROFILE_WINDOW=32`; `just profile GnomonicStars` performs the
locked build, flash, capture, and artifact attestation.
