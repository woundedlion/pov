# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingSpin`).
Raw capture: `build/prof/ringspin_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_ringspin_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses shape scan and ring-SDF `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingSpin profile 70 32` |

Image size:

```text
FLASH: code:52,544, data:147,912, headers:8,440
       free for files:1,822,720
RAM1:  variables:315,072, code:32,344, padding:424
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 737–768 root counter
cycles ÷ 600 MHz match the measured wall sum within **3.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `rs_draw_rings` averages
29.57 ms/frame; its worst window is 35.26 ms/frame
(frames 737–768). Peak frame render is
**49.42 ms** (frames 737–768);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 13.08 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 737–768)

```text
frame                         62.79 ms  37.67 Mcyc  100%
  rs_draw_rings               35.26 ms  21.15 Mcyc   56%
    rs_ring_scan              34.66 ms  20.79 Mcyc   55%
      filter_blend             4.75 ms   2.85 Mcyc    8%  x46520.1 0.1us/c
  rs_timeline_step             76.9 us  46.15 kcyc    0%  x1.0 76.9us/c
  canvas_clear                 87.0 us  52.20 kcyc    0%  x1.0 87.0us/c
  canvas_buffer_wait          27.36 ms  16.42 Mcyc   44%  x1.0 27364.8us/c
```

Wall min/avg/max = 44.59/62.79/80.99 ms. `rs_draw_rings`
accounts for 56.2% of this measured frame at 35.26 ms/frame.
The complete render is 35.42 ms; `canvas_buffer_wait` contributes
27.36 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–32)

```text
frame                         59.16 ms  35.50 Mcyc  100%
  rs_draw_rings               14.96 ms   8.98 Mcyc   25%
    rs_ring_scan              14.54 ms   8.72 Mcyc   25%
      filter_blend             2.08 ms   1.25 Mcyc    4%  x20530.7 0.1us/c
  rs_timeline_step             76.9 us  46.16 kcyc    0%  x1.0 76.9us/c
  canvas_clear                 86.1 us  51.69 kcyc    0%  x1.0 86.1us/c
  canvas_buffer_wait          44.04 ms  26.42 Mcyc   74%  x1.0 44040.1us/c
```

Wall min/avg/max = 0.34/59.16/88.30 ms. `rs_draw_rings`
accounts for 25.3% of this measured frame at 14.96 ms/frame.
The complete render is 15.12 ms; `canvas_buffer_wait` contributes
44.04 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 46,520.1 pixels/frame, 4.49× the 10,368-pixel quadrant. `filter_blend` costs 61.3 cycles/blend; `rs_draw_rings` uses 454.7 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.1/f  0.53/1.65/24.76 us  cpu 3.06%
isr_pack          143.8/f  6.23/6.83/9.54 us  cpu 1.57%
isr_dma_submit    143.8/f  0.61/0.94/2.63 us  cpu 0.22%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.84% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `rs_draw_rings` — 47.4% of aggregated root time, 29.57 ms/frame: inclusive measured scope in the live driver.
2. `rs_ring_scan` — 46.5% of aggregated root time, 28.98 ms/frame: inclusive measured scope in the live driver.
3. `filter_blend` — 6.4% of aggregated root time, 3.99 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches shape scan and ring-SDF `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=RingSpin` and
`HS_PROFILE_WINDOW=32`; `just profile RingSpin` performs the
locked build, flash, capture, and artifact attestation.
