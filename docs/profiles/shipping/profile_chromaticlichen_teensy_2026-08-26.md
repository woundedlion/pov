# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ChromaticLichen`).
Raw capture: `build/prof/chromaticlichen_ship.log`, sourced from the isolated sweep tree.
First archived shipping report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh ChromaticLichen profile 70 32` |

Image size:

```text
FLASH: code:71,280, data:151,664, headers:8,480
       free for files:1,800,192
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 1025–1056 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.4 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
36.14 ms/frame; its worst window is 36.18 ms/frame
(frames 1025–1056). Peak frame render is
**42.87 ms** (frames 257–288);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 19.63 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 1025–1056)

```text
frame                         62.44 ms  37.46 Mcyc  100%
  fx_shader_draw              36.18 ms  21.71 Mcyc   58%  x1.0 36182.6us/c
  fx_prepare_frame             4.11 ms   2.47 Mcyc    7%  x1.0 4109.8us/c
  fx_advance                   2.01 ms   1.21 Mcyc    3%  x1.0 2011.2us/c
  fx_timeline_step             60.1 us  36.08 kcyc    0%  x1.0 60.1us/c
  canvas_clear                 88.8 us  53.27 kcyc    0%  x1.0 88.8us/c
  canvas_buffer_wait          19.96 ms  11.98 Mcyc   32%  x1.0 19962.4us/c
```

Wall min/avg/max = 62.29/62.44/62.58 ms. `fx_shader_draw`
accounts for 57.9% of this measured frame at 36.18 ms/frame.
The complete render is 42.48 ms; `canvas_buffer_wait` contributes
19.96 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 2093.9 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.1/f  0.48/1.68/26.24 us  cpu 3.11%
isr_pack          143.8/f  6.24/6.88/12.62 us  cpu 1.58%
isr_dma_submit    143.8/f  0.59/0.94/11.79 us  cpu 0.22%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.91% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 57.9% of aggregated root time, 36.14 ms/frame: inclusive measured scope in the live driver.
2. `fx_prepare_frame` — 6.5% of aggregated root time, 4.04 ms/frame: inclusive measured scope in the live driver.
3. `fx_advance` — 3.2% of aggregated root time, 1.99 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=ChromaticLichen` and
`HS_PROFILE_WINDOW=32`; `just profile ChromaticLichen` performs the
locked build, flash, capture, and artifact attestation.
