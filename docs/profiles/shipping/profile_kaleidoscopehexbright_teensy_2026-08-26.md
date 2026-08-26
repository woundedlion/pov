# KaleidoscopeHexBright on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexBright`).
Raw capture: `build/prof/kaleidoscopehexbright_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_kaleidoscopehexbright_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexBright 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexBright profile 150 32` |

Image size:

```text
FLASH: code:69,440, data:151,832, headers:9,128
       free for files:1,801,216
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2145–2176 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.8 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
26.48 ms/frame; its worst window is 28.16 ms/frame
(frames 2145–2176). Peak frame render is
**35.55 ms** (frames 1825–1856);
spilled **0/2368 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 26.95 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2145–2176)

```text
frame                         62.45 ms  37.47 Mcyc  100%
  fx_shader_draw              28.16 ms  16.90 Mcyc   45%  x1.0 28162.8us/c
  fx_prepare_frame             3.80 ms   2.28 Mcyc    6%  x1.0 3799.9us/c
  fx_advance                   1.99 ms   1.19 Mcyc    3%  x1.0 1990.1us/c
  fx_timeline_step            103.8 us  62.29 kcyc    0%  x1.0 103.8us/c
  canvas_clear                 89.5 us  53.69 kcyc    0%  x1.0 89.5us/c
  canvas_buffer_wait          28.28 ms  16.97 Mcyc   45%  x1.0 28279.2us/c
```

Wall min/avg/max = 61.01/62.45/63.84 ms. `fx_shader_draw`
accounts for 45.1% of this measured frame at 28.16 ms/frame.
The complete render is 34.17 ms; `canvas_buffer_wait` contributes
28.28 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 641–672)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  fx_shader_draw              23.48 ms  14.09 Mcyc   38%  x1.0 23479.8us/c
  fx_prepare_frame             3.79 ms   2.28 Mcyc    6%  x1.0 3792.6us/c
  fx_advance                   2.04 ms   1.22 Mcyc    3%  x1.0 2035.1us/c
  fx_timeline_step            139.7 us  83.85 kcyc    0%  x1.0 139.7us/c
  canvas_clear                 88.9 us  53.35 kcyc    0%  x1.0 88.9us/c
  canvas_buffer_wait          32.86 ms  19.72 Mcyc   53%  x1.0 32860.0us/c
```

Wall min/avg/max = 61.42/62.42/63.29 ms. `fx_shader_draw`
accounts for 37.6% of this measured frame at 23.48 ms/frame.
The complete render is 29.56 ms; `canvas_buffer_wait` contributes
32.86 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `1` | — | 21/21 | — | 28.16 | 35.55 | 16.0 |
| 2 | `0` | — | 19/19 | — | 28.16 | 34.66 | 16.0 |
| 3 | `2` | — | 34/34 | — | 28.06 | 35.12 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1629.8 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1155.2/f  0.48/1.71/26.84 us  cpu 3.17%
isr_pack          143.9/f  6.24/7.15/11.06 us  cpu 1.65%
isr_dma_submit    143.9/f  0.58/0.93/11.28 us  cpu 0.22%
```

- Pack plus submit costs 1.16 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.03% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 42.4% of aggregated root time, 26.48 ms/frame: inclusive measured scope in the live driver.
2. `fx_prepare_frame` — 6.6% of aggregated root time, 4.14 ms/frame: inclusive measured scope in the live driver.
3. `fx_advance` — 3.2% of aggregated root time, 2.00 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopeHexBright` and
`HS_PROFILE_WINDOW=32`; `just profile KaleidoscopeHexBright` performs the
locked build, flash, capture, and artifact attestation.
