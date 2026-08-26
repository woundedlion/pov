# LatticeMelt on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile LatticeMelt`).
Raw capture: `build/prof/latticemelt_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_latticemelt_teensy_2026-08-18.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | LatticeMelt 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture, `-D HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh LatticeMelt profile 110 16 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size:

```text
FLASH: code:70,856, data:151,672, headers:8,896
       free for files:1,800,192
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 465–480 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.8 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
38.23 ms/frame; its worst window is 38.36 ms/frame
(frames 465–480). Peak frame render is
**43.73 ms** (frames 1–16);
spilled **0/1728 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 18.77 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 465–480)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  fx_shader_draw              38.36 ms  23.02 Mcyc   61%  x1.0 38362.4us/c
  fx_prepare_frame             1.15 ms 689.24 kcyc    2%  x1.0 1148.7us/c
  fx_advance                   2.19 ms   1.31 Mcyc    4%  x1.0 2191.1us/c
  fx_timeline_step             46.0 us  27.61 kcyc    0%  x1.0 46.0us/c
  canvas_clear                 88.1 us  52.86 kcyc    0%  x1.0 88.1us/c
  canvas_buffer_wait          20.56 ms  12.34 Mcyc   33%  x1.0 20558.4us/c
```

Wall min/avg/max = 62.25/62.42/62.62 ms. `fx_shader_draw`
accounts for 61.5% of this measured frame at 38.36 ms/frame.
The complete render is 41.86 ms; `canvas_buffer_wait` contributes
20.56 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `0` | — | 38/38 | — | 38.36 | 43.73 | 16.0 |
| 2 | `2` | — | 67/67 | — | 38.33 | 42.08 | 16.0 |
| 3 | `1` | — | 3/3 | — | 38.30 | 42.00 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 2220.0 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1156.5/f  0.48/1.69/27.26 us  cpu 3.13%
isr_pack          143.9/f  6.24/7.00/21.02 us  cpu 1.61%
isr_dma_submit    143.9/f  0.58/0.94/12.12 us  cpu 0.22%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.96% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 61.3% of aggregated root time, 38.23 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.6% of aggregated root time, 2.22 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.6% of aggregated root time, 0.99 ms/frame: inclusive measured scope in the live driver.

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
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=1200`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=LatticeMelt` and
`HS_PROFILE_WINDOW=16`; `just profile LatticeMelt` performs the
locked build, flash, capture, and artifact attestation.
