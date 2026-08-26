# KaleidoscopePentBright on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopePentBright`).
Raw capture: `build/prof/kaleidoscopepentbright_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_kaleidoscopepentbright_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopePentBright 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopePentBright profile 70 32` |

Image size:

```text
FLASH: code:69,312, data:152,028, headers:9,060
       free for files:1,801,216
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 417–448 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.9 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
28.29 ms/frame; its worst window is 29.94 ms/frame
(frames 417–448). Peak frame render is
**33.37 ms** (frames 865–896);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 29.13 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 417–448)

```text
frame                         62.43 ms  37.46 Mcyc  100%
  fx_shader_draw              29.94 ms  17.96 Mcyc   48%  x1.0 29940.5us/c
  fx_prepare_frame            809.1 us 485.49 kcyc    1%  x1.0 809.1us/c
  fx_advance                   1.98 ms   1.19 Mcyc    3%  x1.0 1977.4us/c
  fx_timeline_step             54.0 us  32.38 kcyc    0%  x1.0 54.0us/c
  canvas_clear                 88.8 us  53.27 kcyc    0%  x1.0 88.8us/c
  canvas_buffer_wait          29.54 ms  17.72 Mcyc   47%  x1.0 29537.8us/c
```

Wall min/avg/max = 61.59/62.43/63.20 ms. `fx_shader_draw`
accounts for 48.0% of this measured frame at 29.94 ms/frame.
The complete render is 32.89 ms; `canvas_buffer_wait` contributes
29.54 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 641–672)

```text
frame                         62.43 ms  37.46 Mcyc  100%
  fx_shader_draw              25.97 ms  15.58 Mcyc   42%  x1.0 25974.3us/c
  fx_prepare_frame            600.1 us 360.04 kcyc    1%  x1.0 600.1us/c
  fx_advance                   1.93 ms   1.16 Mcyc    3%  x1.0 1927.3us/c
  fx_timeline_step             57.6 us  34.54 kcyc    0%  x1.0 57.6us/c
  canvas_clear                 88.1 us  52.85 kcyc    0%  x1.0 88.1us/c
  canvas_buffer_wait          33.76 ms  20.25 Mcyc   54%  x1.0 33756.9us/c
```

Wall min/avg/max = 61.50/62.43/63.48 ms. `fx_shader_draw`
accounts for 41.6% of this measured frame at 25.97 ms/frame.
The complete render is 28.67 ms; `canvas_buffer_wait` contributes
33.76 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1732.7 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1159.0/f  0.48/1.99/23.95 us  cpu 3.69%
isr_pack          143.8/f  6.24/7.07/13.66 us  cpu 1.63%
isr_dma_submit    143.8/f  0.59/0.93/11.18 us  cpu 0.21%
```

- Pack plus submit costs 1.15 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.53% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 45.3% of aggregated root time, 28.29 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.1% of aggregated root time, 1.96 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.4% of aggregated root time, 0.89 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopePentBright` and
`HS_PROFILE_WINDOW=32`; `just profile KaleidoscopePentBright` performs the
locked build, flash, capture, and artifact attestation.
