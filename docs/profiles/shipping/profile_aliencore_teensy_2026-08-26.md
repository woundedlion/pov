# AlienCore on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienCore`).
Raw capture: `build/prof/aliencore_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_aliencore_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienCore 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AlienCore profile 70 32` |

Image size:

```text
FLASH: code:69,304, data:151,664, headers:8,408
       free for files:1,802,240
RAM1:  variables:315,072, code:26,168, padding:6,600
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 33–64 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.3 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
15.72 ms/frame; its worst window is 15.73 ms/frame
(frames 33–64). Peak frame render is
**21.09 ms** (frames 1–32);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 41.41 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 33–64)

```text
frame                         62.44 ms  37.47 Mcyc  100%
  fx_shader_draw              15.73 ms   9.44 Mcyc   25%  x1.0 15733.4us/c
  fx_prepare_frame              1.4 us     856 cyc    0%  x1.0 1.4us/c
  fx_advance                   2.29 ms   1.37 Mcyc    4%  x1.0 2291.2us/c
  fx_timeline_step             59.4 us  35.66 kcyc    0%  x1.0 59.4us/c
  canvas_clear                 88.2 us  52.92 kcyc    0%  x1.0 88.2us/c
  canvas_buffer_wait          44.25 ms  26.55 Mcyc   71%  x1.0 44248.6us/c
```

Wall min/avg/max = 62.30/62.44/62.61 ms. `fx_shader_draw`
accounts for 25.2% of this measured frame at 15.73 ms/frame.
The complete render is 18.19 ms; `canvas_buffer_wait` contributes
44.25 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 910.5 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.7/f  0.48/1.68/30.51 us  cpu 3.12%
isr_pack          143.8/f  6.24/7.10/13.56 us  cpu 1.63%
isr_dma_submit    143.8/f  0.58/0.94/11.63 us  cpu 0.22%
```

- Pack plus submit costs 1.16 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.97% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 25.2% of aggregated root time, 15.72 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.6% of aggregated root time, 2.25 ms/frame: inclusive measured scope in the live driver.
3. `canvas_clear` — 0.1% of aggregated root time, 0.09 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=AlienCore` and
`HS_PROFILE_WINDOW=32`; `just profile AlienCore` performs the
locked build, flash, capture, and artifact attestation.
