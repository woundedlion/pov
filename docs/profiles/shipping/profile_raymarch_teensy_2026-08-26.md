# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Raymarch`).
Raw capture: `build/prof/raymarch_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_raymarch_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses volume scan and SDF-volume `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Raymarch profile 70 32` |

Image size:

```text
FLASH: code:107,904, data:182,948, headers:9,180
       free for files:1,731,584
RAM1:  variables:315,040, code:37,416, padding:28,120
       free for local variables:143,712
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 321–352 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `rm_shader_draw` averages
50.11 ms/frame; its worst window is 52.95 ms/frame
(frames 321–352). Peak frame render is
**62.22 ms** (frames 321–352);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 0.28 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 321–352)

```text
frame                         62.39 ms  37.43 Mcyc  100%
  rm_shader_draw              52.95 ms  31.77 Mcyc   85%
    filter_blend              520.2 us 312.11 kcyc    1%  x5810.1 0.1us/c
  rm_timeline_step            372.2 us 223.30 kcyc    1%  x1.0 372.2us/c
  canvas_clear                 87.3 us  52.38 kcyc    0%  x1.0 87.3us/c
  canvas_buffer_wait           5.86 ms   3.52 Mcyc    9%  x1.0 5864.9us/c
```

Wall min/avg/max = 53.83/62.38/72.37 ms. `rm_shader_draw`
accounts for 84.9% of this measured frame at 52.95 ms/frame.
The complete render is 56.52 ms; `canvas_buffer_wait` contributes
5.86 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 225–256)

```text
frame                         62.46 ms  37.48 Mcyc  100%
  rm_shader_draw              47.15 ms  28.29 Mcyc   75%
    filter_blend              453.2 us 271.91 kcyc    1%  x5008.6 0.1us/c
  rm_timeline_step            377.2 us 226.34 kcyc    1%  x1.0 377.2us/c
  canvas_clear                 86.8 us  52.09 kcyc    0%  x1.0 86.8us/c
  canvas_buffer_wait          11.73 ms   7.04 Mcyc   19%  x1.0 11734.3us/c
```

Wall min/avg/max = 60.66/62.46/64.38 ms. `rm_shader_draw`
accounts for 75.5% of this measured frame at 47.15 ms/frame.
The complete render is 50.73 ms; `canvas_buffer_wait` contributes
11.73 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 5,810.1 pixels/frame, 0.56× the 10,368-pixel quadrant. `filter_blend` costs 53.7 cycles/blend; `rm_shader_draw` uses 5468.2 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.3/f  0.50/1.68/19.43 us  cpu 3.11%
isr_pack          143.9/f  6.26/6.84/9.55 us  cpu 1.57%
isr_dma_submit    143.9/f  0.59/0.93/2.53 us  cpu 0.21%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.90% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `rm_shader_draw` — 80.3% of aggregated root time, 50.11 ms/frame: inclusive measured scope in the live driver.
2. `filter_blend` — 0.8% of aggregated root time, 0.49 ms/frame: inclusive measured scope in the live driver.
3. `rm_timeline_step` — 0.6% of aggregated root time, 0.37 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches volume scan and SDF-volume `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=Raymarch` and
`HS_PROFILE_WINDOW=32`; `just profile Raymarch` performs the
locked build, flash, capture, and artifact attestation.
