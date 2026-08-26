# AlienOcean on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienOcean`).
Raw capture: `build/prof/alienocean_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_alienocean_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienOcean 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienOcean profile 70 32` |

Image size:

```text
FLASH: code:69,384, data:151,664, headers:8,328
       free for files:1,802,240
RAM1:  variables:315,072, code:26,168, padding:6,600
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 97–128 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.9 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
22.18 ms/frame; its worst window is 22.27 ms/frame
(frames 97–128). Peak frame render is
**27.79 ms** (frames 1–32);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 34.71 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 97–128)

```text
frame                         62.45 ms  37.47 Mcyc  100%
  fx_shader_draw              22.27 ms  13.36 Mcyc   36%  x1.0 22273.7us/c
  fx_prepare_frame            552.2 us 331.30 kcyc    1%  x1.0 552.2us/c
  fx_advance                   2.22 ms   1.33 Mcyc    4%  x1.0 2222.0us/c
  fx_timeline_step             51.8 us  31.05 kcyc    0%  x1.0 51.8us/c
  canvas_clear                 88.9 us  53.34 kcyc    0%  x1.0 88.9us/c
  canvas_buffer_wait          37.24 ms  22.34 Mcyc   60%  x1.0 37239.4us/c
```

Wall min/avg/max = 62.27/62.45/62.67 ms. `fx_shader_draw`
accounts for 35.7% of this measured frame at 22.27 ms/frame.
The complete render is 25.21 ms; `canvas_buffer_wait` contributes
37.24 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1289.0 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.9/f  0.48/1.68/26.69 us  cpu 3.11%
isr_pack          143.8/f  6.24/7.00/10.97 us  cpu 1.61%
isr_dma_submit    143.8/f  0.58/0.94/11.16 us  cpu 0.22%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.94% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 35.6% of aggregated root time, 22.18 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.5% of aggregated root time, 2.21 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.0% of aggregated root time, 0.61 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=AlienOcean` and
`HS_PROFILE_WINDOW=32`; `just profile AlienOcean` performs the
locked build, flash, capture, and artifact attestation.
