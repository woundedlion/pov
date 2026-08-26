# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MermaidSkin`).
Raw capture: `build/prof/mermaidskin_ship.log`, sourced from the isolated sweep tree.
First archived shipping report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MermaidSkin 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh MermaidSkin profile 70 32` |

Image size:

```text
FLASH: code:71,416, data:151,660, headers:8,348
       free for files:1,800,192
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 673–704 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.2 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
38.59 ms/frame; its worst window is 38.66 ms/frame
(frames 673–704). Peak frame render is
**45.60 ms** (frames 1–32);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 16.90 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 673–704)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  fx_shader_draw              38.66 ms  23.19 Mcyc   62%  x1.0 38657.8us/c
  fx_prepare_frame             4.03 ms   2.42 Mcyc    6%  x1.0 4027.7us/c
  fx_advance                   1.93 ms   1.16 Mcyc    3%  x1.0 1930.3us/c
  fx_timeline_step             60.2 us  36.16 kcyc    0%  x1.0 60.2us/c
  canvas_clear                 89.0 us  53.42 kcyc    0%  x1.0 89.0us/c
  canvas_buffer_wait          17.64 ms  10.58 Mcyc   28%  x1.0 17638.9us/c
```

Wall min/avg/max = 62.23/62.42/62.60 ms. `fx_shader_draw`
accounts for 61.9% of this measured frame at 38.66 ms/frame.
The complete render is 44.79 ms; `canvas_buffer_wait` contributes
17.64 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 2237.1 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.2/f  0.48/1.74/24.35 us  cpu 3.23%
isr_pack          143.8/f  6.33/7.41/10.01 us  cpu 1.71%
isr_dma_submit    143.8/f  0.59/0.93/11.48 us  cpu 0.22%
```

- Pack plus submit costs 1.20 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.15% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 61.9% of aggregated root time, 38.59 ms/frame: inclusive measured scope in the live driver.
2. `fx_prepare_frame` — 6.9% of aggregated root time, 4.29 ms/frame: inclusive measured scope in the live driver.
3. `fx_advance` — 3.2% of aggregated root time, 1.98 ms/frame: inclusive measured scope in the live driver.

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
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=MermaidSkin` and
`HS_PROFILE_WINDOW=32`; `just profile MermaidSkin` performs the
locked build, flash, capture, and artifact attestation.
