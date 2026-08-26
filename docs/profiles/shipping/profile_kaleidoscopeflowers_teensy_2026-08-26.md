# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeFlowers`).
Raw capture: `build/prof/kaleidoscopeflowers_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_kaleidoscopeflowers_teensy_2026-08-17.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeFlowers profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size:

```text
FLASH: code:70,176, data:151,884, headers:8,340
       free for files:1,801,216
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 1617–1632 root counter
cycles ÷ 600 MHz match the measured wall sum within **4.3 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
30.68 ms/frame; its worst window is 32.94 ms/frame
(frames 1617–1632). Peak frame render is
**36.69 ms** (frames 1601–1616);
spilled **0/4128 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 25.81 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 1617–1632)

```text
frame                         62.43 ms  37.46 Mcyc  100%
  fx_shader_draw              32.94 ms  19.76 Mcyc   53%  x1.0 32939.2us/c
  fx_prepare_frame             1.12 ms 670.72 kcyc    2%  x1.0 1117.8us/c
  fx_advance                   1.94 ms   1.17 Mcyc    3%  x1.0 1944.9us/c
  fx_timeline_step             64.5 us  38.72 kcyc    0%  x1.0 64.5us/c
  canvas_clear                 88.4 us  53.04 kcyc    0%  x1.0 88.4us/c
  canvas_buffer_wait          26.26 ms  15.76 Mcyc   42%  x1.0 26258.7us/c
```

Wall min/avg/max = 61.78/62.43/63.05 ms. `fx_shader_draw`
accounts for 52.8% of this measured frame at 32.94 ms/frame.
The complete render is 36.17 ms; `canvas_buffer_wait` contributes
26.26 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 705–720)

```text
frame                         62.38 ms  37.43 Mcyc  100%
  fx_shader_draw              28.09 ms  16.85 Mcyc   45%  x1.0 28088.2us/c
  fx_prepare_frame            601.3 us 360.82 kcyc    1%  x1.0 601.3us/c
  fx_advance                   1.93 ms   1.16 Mcyc    3%  x1.0 1925.9us/c
  fx_timeline_step            148.4 us  89.09 kcyc    0%  x1.0 148.4us/c
  canvas_clear                 87.6 us  52.60 kcyc    0%  x1.0 87.6us/c
  canvas_buffer_wait          31.51 ms  18.91 Mcyc   51%  x1.0 31512.3us/c
```

Wall min/avg/max = 61.52/62.38/63.15 ms. `fx_shader_draw`
accounts for 45.0% of this measured frame at 28.09 ms/frame.
The complete render is 30.87 ms; `canvas_buffer_wait` contributes
31.51 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `2` | — | 85/85 | — | 32.94 | 36.69 | 16.0 |
| 2 | `1` | — | 67/67 | — | 32.88 | 36.50 | 16.0 |
| 3 | `0` | — | 38/38 | — | 32.79 | 35.92 | 16.0 |
| 4 | `3` | — | 68/68 | — | 32.47 | 36.24 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1906.2 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1153.9/f  0.48/1.71/20.35 us  cpu 3.16%
isr_pack          144.0/f  6.24/7.20/10.01 us  cpu 1.66%
isr_dma_submit    144.0/f  0.58/0.93/12.39 us  cpu 0.22%
```

- Pack plus submit costs 1.17 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.03% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 49.2% of aggregated root time, 30.68 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.2% of aggregated root time, 1.98 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.5% of aggregated root time, 0.94 ms/frame: inclusive measured scope in the live driver.

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
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=2400`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopeFlowers` and
`HS_PROFILE_WINDOW=16`; `just profile KaleidoscopeFlowers` performs the
locked build, flash, capture, and artifact attestation.
