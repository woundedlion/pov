# RingShower on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingShower`).
Raw capture: `build/prof/ringshower_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_ringshower_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses shape scan and ring-SDF `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingShower 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingShower profile 70 32` |

Image size:

```text
FLASH: code:61,688, data:145,636, headers:8,740
       free for files:1,815,552
RAM1:  variables:315,072, code:35,704, padding:29,832
       free for local variables:143,680
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 385–416 root counter
cycles ÷ 600 MHz match the measured wall sum within **0.3 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `rsh_draw_rings` averages
0.85 ms/frame; its worst window is 1.99 ms/frame
(frames 385–416). Peak frame render is
**3.98 ms** (frames 193–224);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 58.52 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 385–416)

```text
frame                         62.39 ms  37.43 Mcyc  100%
  rsh_draw_rings               1.99 ms   1.19 Mcyc    3%
    rsh_ring_plot              1.98 ms   1.19 Mcyc    3%
      filter_blend            165.8 us  99.50 kcyc    0%  x1449.4 0.1us/c
  rsh_timeline_step            54.9 us  32.95 kcyc    0%  x1.0 54.9us/c
  canvas_clear                 84.9 us  50.97 kcyc    0%  x1.0 84.9us/c
  canvas_buffer_wait          60.26 ms  36.16 Mcyc   97%  x1.0 60260.7us/c
```

Wall min/avg/max = 60.33/62.39/63.83 ms. `rsh_draw_rings`
accounts for 3.2% of this measured frame at 1.99 ms/frame.
The complete render is 2.13 ms; `canvas_buffer_wait` contributes
60.26 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–32)

```text
frame                         58.56 ms  35.14 Mcyc  100%
  rsh_draw_rings                0.8 us     482 cyc    0%  x1.0 0.8us/c
  rsh_timeline_step             5.7 us   3.43 kcyc    0%  x1.0 5.7us/c
  canvas_clear                 85.9 us  51.54 kcyc    0%  x1.0 85.9us/c
  canvas_buffer_wait          58.47 ms  35.08 Mcyc  100%  x1.0 58470.5us/c
```

Wall min/avg/max = 0.11/58.56/62.51 ms. `rsh_draw_rings`
accounts for 0.0% of this measured frame at 0.00 ms/frame.
The complete render is 0.09 ms; `canvas_buffer_wait` contributes
58.47 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 1,449.4 pixels/frame, 0.14× the 10,368-pixel quadrant. `filter_blend` costs 68.7 cycles/blend; `rsh_draw_rings` uses 822.1 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.5/f  0.45/1.54/20.18 us  cpu 2.85%
isr_pack          143.7/f  6.23/6.64/9.42 us  cpu 1.53%
isr_dma_submit    143.7/f  0.61/0.94/2.03 us  cpu 0.22%
```

- Pack plus submit costs 1.09 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.60% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `rsh_draw_rings` — 1.4% of aggregated root time, 0.85 ms/frame: inclusive measured scope in the live driver.
2. `rsh_ring_plot` — 1.4% of aggregated root time, 0.85 ms/frame: inclusive measured scope in the live driver.
3. `canvas_clear` — 0.1% of aggregated root time, 0.08 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=RingShower` and
`HS_PROFILE_WINDOW=32`; `just profile RingShower` performs the
locked build, flash, capture, and artifact attestation.
