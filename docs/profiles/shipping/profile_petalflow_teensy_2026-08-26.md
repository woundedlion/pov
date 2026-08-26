# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile PetalFlow`).
Raw capture: `build/prof/petalflow_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_petalflow_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses shape scan and ring-SDF `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh PetalFlow profile 70 32` |

Image size:

```text
FLASH: code:55,824, data:146,972, headers:9,172
       free for files:1,819,648
RAM1:  variables:315,104, code:35,112, padding:30,424
       free for local variables:143,648
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 961–992 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `pf_draw_rings` averages
10.74 ms/frame; its worst window is 11.13 ms/frame
(frames 961–992). Peak frame render is
**11.85 ms** (frames 289–320);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 50.66 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 961–992)

```text
frame                         62.43 ms  37.46 Mcyc  100%
  pf_draw_rings               11.13 ms   6.68 Mcyc   18%
    pf_ring_scan              10.02 ms   6.01 Mcyc   16%
      filter_blend            755.7 us 453.41 kcyc    1%  x6502.7 0.1us/c
    pf_ring_build             998.2 us 598.90 kcyc    2%  x23.2 42.9us/c
  pf_timeline_step             18.4 us  11.04 kcyc    0%  x1.0 18.4us/c
  canvas_clear                 85.0 us  51.03 kcyc    0%  x1.0 85.0us/c
  canvas_buffer_wait          51.19 ms  30.72 Mcyc   82%  x1.0 51193.2us/c
```

Wall min/avg/max = 62.02/62.43/62.89 ms. `pf_draw_rings`
accounts for 17.8% of this measured frame at 11.13 ms/frame.
The complete render is 11.23 ms; `canvas_buffer_wait` contributes
51.19 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 129–160)

```text
frame                         62.46 ms  37.47 Mcyc  100%
  pf_draw_rings               10.08 ms   6.05 Mcyc   16%
    pf_ring_scan               9.05 ms   5.43 Mcyc   14%
      filter_blend            658.4 us 395.05 kcyc    1%  x5732.3 0.1us/c
    pf_ring_build             954.4 us 572.68 kcyc    2%  x22.3 42.8us/c
  pf_timeline_step             18.1 us  10.89 kcyc    0%  x1.0 18.1us/c
  canvas_clear                 85.0 us  50.98 kcyc    0%  x1.0 85.0us/c
  canvas_buffer_wait          52.27 ms  31.36 Mcyc   84%  x1.0 52269.6us/c
```

Wall min/avg/max = 61.67/62.46/63.19 ms. `pf_draw_rings`
accounts for 16.1% of this measured frame at 10.08 ms/frame.
The complete render is 10.19 ms; `canvas_buffer_wait` contributes
52.27 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 6,502.7 pixels/frame, 0.63× the 10,368-pixel quadrant. `filter_blend` costs 69.7 cycles/blend; `pf_draw_rings` uses 1026.9 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.6/f  0.42/1.62/24.64 us  cpu 3.00%
isr_pack          143.8/f  6.24/6.66/9.65 us  cpu 1.53%
isr_dma_submit    143.8/f  0.59/0.94/11.49 us  cpu 0.22%
```

- Pack plus submit costs 1.09 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.74% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `pf_draw_rings` — 17.2% of aggregated root time, 10.74 ms/frame: inclusive measured scope in the live driver.
2. `pf_ring_scan` — 15.5% of aggregated root time, 9.66 ms/frame: inclusive measured scope in the live driver.
3. `pf_ring_build` — 1.6% of aggregated root time, 0.99 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=PetalFlow` and
`HS_PROFILE_WINDOW=32`; `just profile PetalFlow` performs the
locked build, flash, capture, and artifact attestation.
