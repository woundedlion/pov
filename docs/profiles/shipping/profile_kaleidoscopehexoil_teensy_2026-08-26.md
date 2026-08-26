# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexOil`).
Raw capture: `build/prof/kaleidoscopehexoil_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_kaleidoscopehexoil_teensy_2026-08-17.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 140 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexOil profile 140 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size:

```text
FLASH: code:70,208, data:155,772, headers:8,516
       free for files:1,797,120
RAM1:  variables:315,072, code:25,976, padding:6,792
       free for local variables:176,448
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 481–496 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.9 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` averages
32.81 ms/frame; its worst window is 35.18 ms/frame
(frames 481–496). Peak frame render is
**39.04 ms** (frames 513–528);
spilled **0/2208 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 23.46 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 481–496)

```text
frame                         62.43 ms  37.46 Mcyc  100%
  fx_shader_draw              35.18 ms  21.11 Mcyc   56%  x1.0 35178.5us/c
  fx_prepare_frame             1.01 ms 606.89 kcyc    2%  x1.0 1011.4us/c
  fx_advance                   2.27 ms   1.36 Mcyc    4%  x1.0 2274.6us/c
  fx_timeline_step             70.2 us  42.17 kcyc    0%  x1.0 70.2us/c
  canvas_clear                 90.1 us  54.06 kcyc    0%  x1.0 90.1us/c
  canvas_buffer_wait          23.77 ms  14.26 Mcyc   38%  x1.0 23772.6us/c
```

Wall min/avg/max = 62.01/62.43/62.78 ms. `fx_shader_draw`
accounts for 56.4% of this measured frame at 35.18 ms/frame.
The complete render is 38.66 ms; `canvas_buffer_wait` contributes
23.77 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 657–672)

```text
frame                         62.39 ms  37.44 Mcyc  100%
  fx_shader_draw              30.63 ms  18.38 Mcyc   49%  x1.0 30627.8us/c
  fx_prepare_frame            993.9 us 596.36 kcyc    2%  x1.0 993.9us/c
  fx_advance                   2.26 ms   1.36 Mcyc    4%  x1.0 2260.2us/c
  fx_timeline_step            139.6 us  83.75 kcyc    0%  x1.0 139.6us/c
  canvas_clear                 89.8 us  53.85 kcyc    0%  x1.0 89.8us/c
  canvas_buffer_wait          28.25 ms  16.95 Mcyc   45%  x1.0 28248.4us/c
```

Wall min/avg/max = 61.84/62.39/62.72 ms. `fx_shader_draw`
accounts for 49.1% of this measured frame at 30.63 ms/frame.
The complete render is 34.14 ms; `canvas_buffer_wait` contributes
28.25 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `fx_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `fx_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `0` | — | 38/38 | — | 35.18 | 39.04 | 16.0 |
| 2 | `2` | — | 67/67 | — | 35.03 | 38.63 | 16.0 |
| 3 | `1` | — | 33/33 | — | 34.39 | 38.90 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `fx_shader_draw` uses 2035.8 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1155.5/f  0.54/1.70/24.45 us  cpu 3.14%
isr_pack          143.9/f  6.23/6.99/16.00 us  cpu 1.61%
isr_dma_submit    143.9/f  0.59/0.94/11.85 us  cpu 0.22%
```

- Pack plus submit costs 1.14 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.96% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fx_shader_draw` — 52.6% of aggregated root time, 32.81 ms/frame: inclusive measured scope in the live driver.
2. `fx_advance` — 3.6% of aggregated root time, 2.26 ms/frame: inclusive measured scope in the live driver.
3. `fx_prepare_frame` — 1.4% of aggregated root time, 0.86 ms/frame: inclusive measured scope in the live driver.

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

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=KaleidoscopeHexOil` and
`HS_PROFILE_WINDOW=16`; `just profile KaleidoscopeHexOil` performs the
locked build, flash, capture, and artifact attestation.
