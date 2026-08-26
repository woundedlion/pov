# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeSmooth`).
Raw capture: `build/prof/kaleidoscopesmooth_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_kaleidoscopesmooth_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeSmooth profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size:

```text
FLASH: code:69,968, data:151,884, headers:8,548
       free for files:1,801,216
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 3889–3904 root counter
cycles ÷ 600 MHz match the measured wall sum within **0.5 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
27.09 ms/frame; its worst window is 32.30 ms/frame
(frames 3889–3904). Peak frame render is
**35.85 ms** (frames 1041–1056);
spilled **0/4128 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 26.65 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 3889–3904)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  fx_shader_draw              32.30 ms  19.38 Mcyc   52%  x1.0 32296.0us/c
  fx_prepare_frame            631.6 us 378.98 kcyc    1%  x1.0 631.6us/c
  fx_advance                   2.00 ms   1.20 Mcyc    3%  x1.0 2004.4us/c
  fx_timeline_step            149.8 us  89.88 kcyc    0%  x1.0 149.8us/c
  canvas_clear                 87.8 us  52.69 kcyc    0%  x1.0 87.8us/c
  canvas_buffer_wait          27.24 ms  16.34 Mcyc   44%  x1.0 27237.6us/c
```

Wall min/avg/max = 61.30/62.42/63.63 ms. `fx_shader_draw`
accounts for 51.7% of this measured frame at 32.30 ms/frame.
The complete render is 35.19 ms; `canvas_buffer_wait` contributes
27.24 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 2481–2496)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  fx_shader_draw              24.04 ms  14.42 Mcyc   39%  x1.0 24037.4us/c
  fx_prepare_frame            824.6 us 494.81 kcyc    1%  x1.0 824.6us/c
  fx_advance                   2.06 ms   1.23 Mcyc    3%  x1.0 2056.4us/c
  fx_timeline_step             76.1 us  45.64 kcyc    0%  x1.0 76.1us/c
  canvas_clear                 87.8 us  52.67 kcyc    0%  x1.0 87.8us/c
  canvas_buffer_wait          35.33 ms  21.20 Mcyc   57%  x1.0 35325.1us/c
```

Wall min/avg/max = 61.96/62.42/62.87 ms. `fx_shader_draw`
accounts for 38.5% of this measured frame at 24.04 ms/frame.
The complete render is 27.10 ms; `canvas_buffer_wait` contributes
35.33 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `1` | — | 18/18 | — | 32.30 | 35.81 | 16.0 |
| 2 | `2` | — | 67/67 | — | 31.97 | 35.85 | 16.0 |
| 3 | `0` | — | 38/38 | — | 30.73 | 34.20 | 16.0 |
| 4 | `4` | — | 67/67 | — | 27.92 | 31.50 | 16.0 |
| 5 | `3` | — | 68/68 | — | 27.37 | 31.33 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1869.0 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1153.9/f  0.48/1.69/24.55 us  cpu 3.13%
isr_pack          144.0/f  6.24/6.99/10.03 us  cpu 1.61%
isr_dma_submit    144.0/f  0.58/0.93/11.94 us  cpu 0.22%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.95% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 43.4% of aggregated root time, 27.09 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.2% of aggregated root time, 2.01 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.5% of aggregated root time, 0.95 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopeSmooth` and
`HS_PROFILE_WINDOW=16`; `just profile KaleidoscopeSmooth` performs the
locked build, flash, capture, and artifact attestation.
