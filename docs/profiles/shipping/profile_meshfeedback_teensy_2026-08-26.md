# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `build/prof/meshfeedback_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_meshfeedback_teensy_2026-08-05.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses pixel-feedback, filter-pipeline, screen-AA, and scan `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture, `-D HS_PROFILE_EPOCH_REVS=3400` |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"` |

Image size:

```text
FLASH: code:116,800, data:186,712, headers:8,808
       free for files:1,719,296
RAM1:  variables:315,232, code:63,752, padding:1,784
       free for local variables:143,520
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 3921–3936 root counter
cycles ÷ 600 MHz match the measured wall sum within **0.3 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `mf_feedback_flush` averages
41.37 ms/frame; its worst window is 47.23 ms/frame
(frames 3921–3936). Peak frame render is
**58.30 ms** (frames 1201–1216);
spilled **0/6688 frames**
(0.0%).

A display window is 62.5 ms. This effect renders the **FULL canvas**, 288×144 = 41,472 pixels, four times a quadrant. The peak retains 4.20 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 3921–3936)

```text
frame                         62.27 ms  37.36 Mcyc  100%
  mf_mesh_draw                 5.71 ms   3.43 Mcyc    9%
    filter_blend               1.02 ms 611.96 kcyc    2%  x9227.1 0.1us/c
  mf_feedback_flush           47.23 ms  28.34 Mcyc   76%
    feedback_composite        39.93 ms  23.96 Mcyc   64%  x1.0 39930.0us/c
    feedback_populate          7.27 ms   4.36 Mcyc   12%  x1.0 7266.0us/c
    feedback_litscan            0.1 us     109 cyc    0%  x1.0 0.1us/c
  mf_timeline_step             49.8 us  29.91 kcyc    0%  x1.0 49.8us/c
  mf_apply_params               0.2 us     118 cyc    0%  x1.0 0.2us/c
  canvas_clear                408.9 us 245.39 kcyc    1%  x1.0 408.9us/c
  canvas_buffer_wait           8.87 ms   5.32 Mcyc   14%  x1.0 8869.4us/c
```

Wall min/avg/max = 60.24/62.27/63.16 ms. `mf_feedback_flush`
accounts for 75.8% of this measured frame at 47.23 ms/frame.
The complete render is 53.40 ms; `canvas_buffer_wait` contributes
8.87 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 5585–5600)

```text
frame                         62.44 ms  37.46 Mcyc  100%
  mf_mesh_draw                 9.14 ms   5.48 Mcyc   15%
    filter_blend               1.59 ms 951.43 kcyc    3%  x14445.8 0.1us/c
  mf_feedback_flush           33.03 ms  19.82 Mcyc   53%
    feedback_composite        33.01 ms  19.80 Mcyc   53%  x1.0 33007.9us/c
    feedback_populate           0.0 us      35 cyc    0%  x1.0 0.0us/c
    feedback_litscan            0.1 us     111 cyc    0%  x1.0 0.1us/c
  mf_timeline_step             46.6 us  27.97 kcyc    0%  x1.0 46.6us/c
  mf_apply_params               0.2 us     114 cyc    0%  x1.0 0.2us/c
  canvas_clear                406.8 us 244.07 kcyc    1%  x1.0 406.8us/c
  canvas_buffer_wait          19.82 ms  11.89 Mcyc   32%  x1.0 19815.8us/c
```

Wall min/avg/max = 61.34/62.44/63.14 ms. `mf_feedback_flush`
accounts for 52.9% of this measured frame at 33.03 ms/frame.
The complete render is 42.62 ms; `canvas_buffer_wait` contributes
19.82 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `mf_feedback_flush` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `mf_feedback_flush` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `5` | — | 30/30 | 9,227.1 | 47.23 | 56.71 | 16.0 |
| 2 | `3` | — | 45/45 | 9,610.2 | 45.91 | 55.64 | 16.0 |
| 3 | `4` | — | 41/41 | 3,713.4 | 45.03 | 52.19 | 16.0 |
| 4 | `2` | — | 46/46 | 5,291.2 | 45.01 | 50.70 | 16.0 |
| 5 | `9` | — | 30/30 | 10,565.7 | 44.89 | 55.91 | 16.0 |
| 6 | `8` | — | 30/30 | 10,249.6 | 43.83 | 50.64 | 16.0 |
| 7 | `6` | — | 31/31 | 8,899.2 | 43.35 | 58.30 | 16.0 |
| 8 | `7` | — | 30/30 | 8,265.4 | 42.30 | 49.30 | 16.0 |
| 9 | `10` | — | 30/30 | 8,214.3 | 42.25 | 49.83 | 16.0 |
| 10 | `1` | — | 30/30 | 6,701.2 | 40.91 | 51.60 | 16.0 |
| 11 | `11` | — | 30/30 | 11,905.6 | 40.75 | 55.06 | 16.0 |
| 12 | `0` | — | 15/15 | 8,113.5 | 40.39 | 47.02 | 16.0 |
| 13 | `12` | — | 30/30 | 14,429.9 | 39.02 | 55.93 | 16.0 |

### Per-pixel figures

The peak dominant-scope window blends 9,227.1 pixels/frame, 0.89× the 10,368-pixel quadrant. `filter_blend` costs 66.3 cycles/blend; `mf_feedback_flush` uses 3070.9 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1153.2/f  0.43/1.73/19.62 us  cpu 3.20%
isr_pack          144.0/f  6.24/7.16/10.07 us  cpu 1.65%
isr_dma_submit    144.0/f  0.65/0.95/9.83 us  cpu 0.22%
```

- Pack plus submit costs 1.17 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.06% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `mf_feedback_flush` — 66.3% of aggregated root time, 41.37 ms/frame: inclusive measured scope in the live driver.
2. `feedback_composite` — 56.9% of aggregated root time, 35.50 ms/frame: inclusive measured scope in the live driver.
3. `feedback_populate` — 9.3% of aggregated root time, 5.84 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches pixel-feedback, filter-pipeline, screen-AA, and scan `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=3400`.
- Provenance attests a clean source tree at `20ca3cb48892795a4575dd9d16d31a699f79df75`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=MeshFeedback` and
`HS_PROFILE_WINDOW=16`; `just profile MeshFeedback` performs the
locked build, flash, capture, and artifact attestation.
