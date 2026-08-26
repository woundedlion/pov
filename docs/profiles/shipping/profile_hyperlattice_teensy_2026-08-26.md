# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Raw capture: `build/prof/hyperlattice_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_hyperlattice_teensy_2026-08-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses `Scan::Shader::draw_cached` plus the inlined lens, noise, and OKLab/gamut `HS_O3_FN` helpers used by this fixed pipeline |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size:

```text
FLASH: code:62,232, data:147,420, headers:8,460
       free for files:1,813,504
RAM1:  variables:315,008, code:22,312, padding:10,456
       free for local variables:176,512
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 849–864 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.9 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hl_shader_draw` averages
37.68 ms/frame; its worst window is 46.07 ms/frame
(frames 849–864). Peak frame render is
**51.19 ms** (frames 545–560);
spilled **0/2688 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 11.31 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 849–864)

```text
frame                         62.50 ms  37.50 Mcyc  100%
  hl_shader_draw              46.07 ms  27.64 Mcyc   74%  x1.0 46074.0us/c
  hl_timeline_step              1.9 us   1.19 kcyc    0%  x1.0 1.9us/c
  canvas_clear                 87.4 us  52.48 kcyc    0%  x1.0 87.4us/c
  canvas_buffer_wait          13.84 ms   8.30 Mcyc   22%  x1.0 13838.4us/c
```

Wall min/avg/max = 60.14/62.50/64.98 ms. `hl_shader_draw`
accounts for 73.7% of this measured frame at 46.07 ms/frame.
The complete render is 48.66 ms; `canvas_buffer_wait` contributes
13.84 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1233–1248)

```text
frame                         62.40 ms  37.44 Mcyc  100%
  hl_shader_draw              28.12 ms  16.87 Mcyc   45%  x1.0 28123.4us/c
  hl_timeline_step              1.9 us   1.18 kcyc    0%  x1.0 1.9us/c
  canvas_clear                 85.6 us  51.40 kcyc    0%  x1.0 85.6us/c
  canvas_buffer_wait          31.53 ms  18.92 Mcyc   51%  x1.0 31525.4us/c
```

Wall min/avg/max = 60.79/62.40/64.13 ms. `hl_shader_draw`
accounts for 45.1% of this measured frame at 28.12 ms/frame.
The complete render is 30.87 ms; `canvas_buffer_wait` contributes
31.53 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `hl_shader_draw` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `hl_shader_draw` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `2` | — | 70/70 | — | 46.07 | 51.19 | 16.0 |
| 2 | `3` | — | 43/43 | — | 44.00 | 48.95 | 16.0 |
| 3 | `0` | — | 20/20 | — | 37.81 | 42.16 | 16.0 |
| 4 | `1` | — | 35/35 | — | 37.73 | 41.77 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `hl_shader_draw` uses 2666.3 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.9/f  0.45/1.65/20.48 us  cpu 3.05%
isr_pack          143.9/f  6.23/6.77/16.70 us  cpu 1.56%
isr_dma_submit    143.9/f  0.59/0.95/12.27 us  cpu 0.22%
```

- Pack plus submit costs 1.11 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.83% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `hl_shader_draw` — 60.4% of aggregated root time, 37.68 ms/frame: inclusive measured scope in the live driver.
2. `canvas_clear` — 0.1% of aggregated root time, 0.09 ms/frame: inclusive measured scope in the live driver.
3. `hl_timeline_step` — 0.0% of aggregated root time, 0.01 ms/frame: inclusive measured scope in the live driver.

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
- Provenance attests a clean source tree at `20ca3cb48892795a4575dd9d16d31a699f79df75`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HyperLattice` and
`HS_PROFILE_WINDOW=16`; `just profile HyperLattice` performs the
locked build, flash, capture, and artifact attestation.
