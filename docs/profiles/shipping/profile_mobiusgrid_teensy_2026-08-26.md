# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MobiusGrid`).
Raw capture: `build/prof/mobiusgrid_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_mobiusgrid_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh MobiusGrid profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size:

```text
FLASH: code:70,008, data:152,304, headers:9,112
       free for files:1,800,192
RAM1:  variables:315,072, code:26,232, padding:6,536
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2001–2016 root counter
cycles ÷ 600 MHz match the measured wall sum within **3.0 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
18.67 ms/frame; its worst window is 19.08 ms/frame
(frames 2001–2016). Peak frame render is
**22.37 ms** (frames 785–800);
spilled **0/2688 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 40.13 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2001–2016)

```text
frame                         62.44 ms  37.47 Mcyc  100%
  fx_shader_draw              19.08 ms  11.45 Mcyc   31%  x1.0 19083.1us/c
  fx_prepare_frame            700.4 us 420.28 kcyc    1%  x1.0 700.4us/c
  fx_advance                   2.25 ms   1.35 Mcyc    4%  x1.0 2250.9us/c
  fx_timeline_step            157.4 us  94.47 kcyc    0%  x1.0 157.4us/c
  canvas_clear                 87.8 us  52.66 kcyc    0%  x1.0 87.8us/c
  canvas_buffer_wait          40.14 ms  24.09 Mcyc   64%  x1.0 40144.1us/c
```

Wall min/avg/max = 62.19/62.44/62.58 ms. `fx_shader_draw`
accounts for 30.6% of this measured frame at 19.08 ms/frame.
The complete render is 22.30 ms; `canvas_buffer_wait` contributes
40.14 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–16)

```text
frame                         57.31 ms  34.39 Mcyc  100%
  fx_shader_draw              18.44 ms  11.06 Mcyc   32%  x1.0 18438.4us/c
  fx_prepare_frame            563.9 us 338.33 kcyc    1%  x1.0 563.9us/c
  fx_advance                   2.35 ms   1.41 Mcyc    4%  x1.0 2347.2us/c
  fx_timeline_step             54.1 us  32.47 kcyc    0%  x1.0 54.1us/c
  canvas_clear                 89.8 us  53.92 kcyc    0%  x1.0 89.8us/c
  canvas_buffer_wait          35.79 ms  21.48 Mcyc   62%  x1.0 35794.6us/c
```

Wall min/avg/max = 21.08/57.31/62.80 ms. `fx_shader_draw`
accounts for 32.2% of this measured frame at 18.44 ms/frame.
The complete render is 21.51 ms; `canvas_buffer_wait` contributes
35.79 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `1` | — | 63/63 | — | 19.08 | 22.37 | 16.0 |
| 2 | `2` | — | 67/67 | — | 19.07 | 22.37 | 16.0 |
| 3 | `0` | — | 38/38 | — | 18.78 | 21.73 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1104.3 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.8/f  0.48/1.66/31.84 us  cpu 3.07%
isr_pack          143.9/f  6.25/6.80/17.72 us  cpu 1.57%
isr_dma_submit    143.9/f  0.58/0.94/11.68 us  cpu 0.22%
```

- Pack plus submit costs 1.11 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.85% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 29.9% of aggregated root time, 18.67 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.6% of aggregated root time, 2.23 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.0% of aggregated root time, 0.60 ms/frame: inclusive measured scope in the live driver.

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
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=1600`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=MobiusGrid` and
`HS_PROFILE_WINDOW=16`; `just profile MobiusGrid` performs the
locked build, flash, capture, and artifact attestation.
