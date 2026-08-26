# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AlienBrain`).
Raw capture: `build/prof/alienbrain_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_alienbrain_teensy_2026-08-16.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AlienBrain 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 300 s capture, `-D HS_PROFILE_EPOCH_REVS=2600` |
| Reproduce | `bash tools/profile_one.sh AlienBrain profile 300 16 "-D HS_PROFILE_EPOCH_REVS=2600"` |

Image size:

```text
FLASH: code:69,776, data:151,808, headers:8,816
       free for files:1,801,216
RAM1:  variables:315,072, code:26,136, padding:6,632
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 1441–1456 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
25.44 ms/frame; its worst window is 26.03 ms/frame
(frames 1441–1456). Peak frame render is
**31.61 ms** (frames 1–16);
spilled **0/4768 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 30.89 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 1441–1456)

```text
frame                         62.44 ms  37.47 Mcyc  100%
  fx_shader_draw              26.03 ms  15.62 Mcyc   42%  x1.0 26027.9us/c
  fx_prepare_frame            660.6 us 396.40 kcyc    1%  x1.0 660.6us/c
  fx_advance                   2.22 ms   1.33 Mcyc    4%  x1.0 2217.2us/c
  fx_timeline_step             67.4 us  40.50 kcyc    0%  x1.0 67.4us/c
  canvas_clear                 89.0 us  53.42 kcyc    0%  x1.0 89.0us/c
  canvas_buffer_wait          33.37 ms  20.02 Mcyc   53%  x1.0 33368.8us/c
```

Wall min/avg/max = 62.23/62.44/62.58 ms. `fx_shader_draw`
accounts for 41.7% of this measured frame at 26.03 ms/frame.
The complete render is 29.07 ms; `canvas_buffer_wait` contributes
33.37 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 2465–2480)

```text
frame                         62.44 ms  37.46 Mcyc  100%
  fx_shader_draw              24.38 ms  14.63 Mcyc   39%  x1.0 24384.9us/c
  fx_prepare_frame            858.2 us 514.96 kcyc    1%  x1.0 858.2us/c
  fx_advance                   2.21 ms   1.33 Mcyc    4%  x1.0 2210.6us/c
  fx_timeline_step             68.4 us  41.09 kcyc    0%  x1.0 68.4us/c
  canvas_clear                 89.0 us  53.43 kcyc    0%  x1.0 89.0us/c
  canvas_buffer_wait          34.81 ms  20.89 Mcyc   56%  x1.0 34810.7us/c
```

Wall min/avg/max = 62.27/62.44/62.56 ms. `fx_shader_draw`
accounts for 39.1% of this measured frame at 24.38 ms/frame.
The complete render is 27.62 ms; `canvas_buffer_wait` contributes
34.81 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `2` | — | 67/67 | — | 26.03 | 29.47 | 16.0 |
| 2 | `0` | — | 38/38 | — | 25.97 | 31.61 | 16.0 |
| 3 | `1` | — | 58/58 | — | 25.97 | 29.34 | 16.0 |
| 4 | `4` | — | 67/67 | — | 25.93 | 29.56 | 16.0 |
| 5 | `3` | — | 68/68 | — | 25.26 | 29.12 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 1506.2 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1153.6/f  0.48/1.67/28.45 us  cpu 3.07%
isr_pack          144.0/f  6.24/6.81/19.92 us  cpu 1.57%
isr_dma_submit    144.0/f  0.59/0.94/11.69 us  cpu 0.22%
```

- Pack plus submit costs 1.11 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.86% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 40.8% of aggregated root time, 25.44 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.6% of aggregated root time, 2.24 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.3% of aggregated root time, 0.84 ms/frame: inclusive measured scope in the live driver.

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
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=2600`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=AlienBrain` and
`HS_PROFILE_WINDOW=16`; `just profile AlienBrain` performs the
locked build, flash, capture, and artifact attestation.
