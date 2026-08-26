# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GridSpace`).
Raw capture: `build/prof/gridspace_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_gridspace_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GridSpace 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GridSpace profile 70 32` |

Image size:

```text
FLASH: code:70,056, data:152,220, headers:9,148
       free for files:1,800,192
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 577–608 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.1 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
21.62 ms/frame; its worst window is 21.68 ms/frame
(frames 577–608). Peak frame render is
**28.38 ms** (frames 417–448);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 34.12 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 577–608)

```text
frame                         62.44 ms  37.47 Mcyc  100%
  fx_shader_draw              21.68 ms  13.01 Mcyc   35%  x1.0 21679.0us/c
  fx_prepare_frame             4.05 ms   2.43 Mcyc    6%  x1.0 4052.8us/c
  fx_advance                   2.21 ms   1.33 Mcyc    4%  x1.0 2211.1us/c
  fx_timeline_step             63.7 us  38.24 kcyc    0%  x1.0 63.7us/c
  canvas_clear                 86.7 us  52.03 kcyc    0%  x1.0 86.7us/c
  canvas_buffer_wait          34.29 ms  20.58 Mcyc   55%  x1.0 34293.7us/c
```

Wall min/avg/max = 59.37/62.44/65.55 ms. `fx_shader_draw`
accounts for 34.7% of this measured frame at 21.68 ms/frame.
The complete render is 28.15 ms; `canvas_buffer_wait` contributes
34.29 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1254.6 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.9/f  0.50/1.71/28.82 us  cpu 3.16%
isr_pack          143.8/f  6.24/6.93/15.56 us  cpu 1.60%
isr_dma_submit    143.8/f  0.59/0.94/10.99 us  cpu 0.22%
```

- Pack plus submit costs 1.13 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.97% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 34.7% of aggregated root time, 21.62 ms/frame: inclusive measured scope in the live driver.
2. `fx_prepare_frame` — 6.4% of aggregated root time, 3.98 ms/frame: inclusive measured scope in the live driver.
3. `fx_advance` — 3.6% of aggregated root time, 2.24 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=GridSpace` and
`HS_PROFILE_WINDOW=32`; `just profile GridSpace` performs the
locked build, flash, capture, and artifact attestation.
