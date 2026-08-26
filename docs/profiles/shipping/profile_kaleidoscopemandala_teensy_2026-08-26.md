# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeMandala`).
Raw capture: `build/prof/kaleidoscopemandala_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_kaleidoscopemandala_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeMandala 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeMandala profile 150 32` |

Image size:

```text
FLASH: code:69,872, data:151,964, headers:8,564
       free for files:1,801,216
RAM1:  variables:315,072, code:26,168, padding:6,600
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 1601–1632 root counter
cycles ÷ 600 MHz match the measured wall sum within **3.9 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
30.61 ms/frame; its worst window is 33.08 ms/frame
(frames 1601–1632). Peak frame render is
**40.29 ms** (frames 1953–1984);
spilled **0/2368 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 22.21 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 1601–1632)

```text
frame                         62.44 ms  37.46 Mcyc  100%
  fx_shader_draw              33.08 ms  19.85 Mcyc   53%  x1.0 33083.6us/c
  fx_prepare_frame            930.3 us 558.21 kcyc    1%  x1.0 930.3us/c
  fx_advance                   2.23 ms   1.34 Mcyc    4%  x1.0 2230.1us/c
  fx_timeline_step             51.8 us  31.12 kcyc    0%  x1.0 51.8us/c
  canvas_clear                 90.4 us  54.23 kcyc    0%  x1.0 90.4us/c
  canvas_buffer_wait          26.03 ms  15.62 Mcyc   42%  x1.0 26034.6us/c
```

Wall min/avg/max = 61.36/62.44/63.54 ms. `fx_shader_draw`
accounts for 53.0% of this measured frame at 33.08 ms/frame.
The complete render is 36.41 ms; `canvas_buffer_wait` contributes
26.03 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 289–320)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  fx_shader_draw              28.14 ms  16.88 Mcyc   45%  x1.0 28135.7us/c
  fx_prepare_frame             1.12 ms 670.92 kcyc    2%  x1.0 1118.2us/c
  fx_advance                   2.22 ms   1.33 Mcyc    4%  x1.0 2220.9us/c
  fx_timeline_step             52.7 us  31.62 kcyc    0%  x1.0 52.7us/c
  canvas_clear                 90.1 us  54.04 kcyc    0%  x1.0 90.1us/c
  canvas_buffer_wait          30.78 ms  18.47 Mcyc   49%  x1.0 30777.4us/c
```

Wall min/avg/max = 61.34/62.42/63.62 ms. `fx_shader_draw`
accounts for 45.1% of this measured frame at 28.14 ms/frame.
The complete render is 31.64 ms; `canvas_buffer_wait` contributes
30.78 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `2` | — | 34/34 | — | 33.08 | 40.04 | 16.0 |
| 2 | `1` | — | 21/21 | — | 32.56 | 40.29 | 16.0 |
| 3 | `0` | — | 19/19 | — | 31.89 | 35.91 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1914.6 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1155.3/f  0.48/1.76/29.25 us  cpu 3.25%
isr_pack          143.9/f  6.24/7.26/15.29 us  cpu 1.67%
isr_dma_submit    143.9/f  0.58/0.94/11.15 us  cpu 0.22%
```

- Pack plus submit costs 1.18 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.14% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 49.1% of aggregated root time, 30.61 ms/frame: inclusive measured scope in the live driver.
2. `fx_prepare_frame` — 3.8% of aggregated root time, 2.35 ms/frame: inclusive measured scope in the live driver.
3. `fx_advance` — 3.6% of aggregated root time, 2.25 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `none`.
- Provenance attests a clean source tree at `20ca3cb48892795a4575dd9d16d31a699f79df75`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopeMandala` and
`HS_PROFILE_WINDOW=32`; `just profile KaleidoscopeMandala` performs the
locked build, flash, capture, and artifact attestation.
