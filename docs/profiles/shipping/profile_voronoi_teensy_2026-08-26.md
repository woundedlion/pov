# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Voronoi`).
Raw capture: `build/prof/voronoi_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_voronoi_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses shader scan, shading, and math `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Voronoi 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Voronoi profile 70 32` |

Image size:

```text
FLASH: code:39,344, data:145,824, headers:8,368
       free for files:1,838,080
RAM1:  variables:315,040, code:20,472, padding:12,296
       free for local variables:176,480
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 353–384 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.7 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `vo_shade` averages
8.22 ms/frame; its worst window is 8.30 ms/frame
(frames 353–384). Peak frame render is
**8.96 ms** (frames 993–1024);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 53.54 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 353–384)

```text
frame                         62.46 ms  37.47 Mcyc  100%
  vo_shade                     8.30 ms   4.98 Mcyc   13%  x1.0 8299.9us/c
  vo_kdtree                   342.0 us 205.19 kcyc    1%  x1.0 342.0us/c
  vo_animate                   51.0 us  30.62 kcyc    0%  x1.0 51.0us/c
  canvas_clear                 85.0 us  51.00 kcyc    0%  x1.0 85.0us/c
  canvas_buffer_wait          53.68 ms  32.21 Mcyc   86%  x1.0 53679.1us/c
```

Wall min/avg/max = 62.23/62.46/62.71 ms. `vo_shade`
accounts for 13.3% of this measured frame at 8.30 ms/frame.
The complete render is 8.78 ms; `canvas_buffer_wait` contributes
53.68 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 289–320)

```text
frame                         62.46 ms  37.47 Mcyc  100%
  vo_shade                     8.14 ms   4.88 Mcyc   13%  x1.0 8139.4us/c
  vo_kdtree                   339.8 us 203.87 kcyc    1%  x1.0 339.8us/c
  vo_animate                   51.0 us  30.62 kcyc    0%  x1.0 51.0us/c
  canvas_clear                 84.9 us  50.95 kcyc    0%  x1.0 84.9us/c
  canvas_buffer_wait          53.84 ms  32.30 Mcyc   86%  x1.0 53840.6us/c
```

Wall min/avg/max = 62.17/62.46/62.71 ms. `vo_shade`
accounts for 13.0% of this measured frame at 8.14 ms/frame.
The complete render is 8.62 ms; `canvas_buffer_wait` contributes
53.84 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `vo_shade` uses 480.3 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.6/f  0.45/1.65/17.36 us  cpu 3.05%
isr_pack          143.8/f  6.28/6.96/9.62 us  cpu 1.60%
isr_dma_submit    143.8/f  0.62/0.95/2.16 us  cpu 0.22%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.87% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `vo_shade` — 13.2% of aggregated root time, 8.22 ms/frame: inclusive measured scope in the live driver.
2. `vo_kdtree` — 0.6% of aggregated root time, 0.35 ms/frame: inclusive measured scope in the live driver.
3. `canvas_clear` — 0.1% of aggregated root time, 0.09 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches shader scan, shading, and math `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=Voronoi` and
`HS_PROFILE_WINDOW=32`; `just profile Voronoi` performs the
locked build, flash, capture, and artifact attestation.
