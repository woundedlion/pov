# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `build/prof/shapeshifter_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_shapeshifter_teensy_2026-08-08.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses Plot cull/raster and shape-scan `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size:

```text
FLASH: code:76,976, data:150,696, headers:8,872
       free for files:1,795,072
RAM1:  variables:315,072, code:48,744, padding:16,792
       free for local variables:143,680
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2353–2368 root counter
cycles ÷ 600 MHz match the measured wall sum within **3.1 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `ss_draw_all` averages
23.47 ms/frame; its worst window is 51.53 ms/frame
(frames 2353–2368). Peak frame render is
**58.13 ms** (frames 2353–2368);
spilled **0/2448 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 4.37 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2353–2368)

```text
frame                         62.42 ms  37.45 Mcyc  100%
  ss_draw_all                 51.53 ms  30.92 Mcyc   83%
    ss_plot_dispatch          51.35 ms  30.81 Mcyc   82%  x136.9 375.1us/c
  ss_timeline_step             53.4 us  32.03 kcyc    0%  x1.0 53.4us/c
  ss_buffer_wait              10.83 ms   6.50 Mcyc   17%
    canvas_clear               88.8 us  53.32 kcyc    0%  x1.0 88.8us/c
    canvas_buffer_wait        10.75 ms   6.45 Mcyc   17%  x1.0 10745.1us/c
```

Wall min/avg/max = 50.88/62.42/74.48 ms. `ss_draw_all`
accounts for 82.5% of this measured frame at 51.53 ms/frame.
The complete render is 51.59 ms; `ss_buffer_wait` contributes
10.83 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 545–560)

```text
frame                         62.50 ms  37.50 Mcyc  100%
  ss_draw_all                  5.10 ms   3.06 Mcyc    8%
    ss_plot_dispatch           4.99 ms   2.99 Mcyc    8%  x24.2 206.3us/c
  ss_timeline_step             55.7 us  33.44 kcyc    0%  x1.0 55.7us/c
  ss_buffer_wait              57.34 ms  34.41 Mcyc   92%
    canvas_clear               84.5 us  50.73 kcyc    0%  x1.0 84.5us/c
    canvas_buffer_wait        57.26 ms  34.35 Mcyc   92%  x1.0 57257.4us/c
```

Wall min/avg/max = 61.10/62.50/63.91 ms. `ss_draw_all`
accounts for 8.2% of this measured frame at 5.10 ms/frame.
The complete render is 5.16 ms; `ss_buffer_wait` contributes
57.34 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `ss_draw_all` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `ss_draw_all` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `1` | — | 15/15 | — | 51.53 | 58.13 | 16.0 |
| 2 | `0` | — | 15/15 | — | 45.66 | 52.62 | 16.0 |
| 3 | `9` | — | 15/15 | — | 36.83 | 37.46 | 16.0 |
| 4 | `4` | — | 15/15 | — | 35.53 | 37.13 | 16.0 |
| 5 | `8` | — | 15/15 | — | 24.35 | 26.69 | 16.0 |
| 6 | `7` | — | 15/15 | — | 19.91 | 21.77 | 16.0 |
| 7 | `6` | — | 15/15 | — | 17.39 | 18.05 | 16.0 |
| 8 | `2` | — | 18/18 | — | 12.88 | 15.23 | 16.0 |
| 9 | `5` | — | 15/15 | — | 11.23 | 13.64 | 16.0 |
| 10 | `3` | — | 15/15 | — | 9.75 | 9.86 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `ss_draw_all` uses 2982.0 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1155.0/f  0.45/1.66/20.44 us  cpu 3.06%
isr_pack          143.9/f  6.24/6.87/9.65 us  cpu 1.58%
isr_dma_submit    143.9/f  0.58/0.93/3.02 us  cpu 0.21%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.86% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `ss_draw_all` — 37.6% of aggregated root time, 23.47 ms/frame: inclusive measured scope in the live driver.
2. `ss_plot_dispatch` — 37.4% of aggregated root time, 23.36 ms/frame: inclusive measured scope in the live driver.
3. `canvas_clear` — 0.1% of aggregated root time, 0.09 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches Plot cull/raster and shape-scan `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=1600`.
- Provenance attests a clean source tree at `20ca3cb48892795a4575dd9d16d31a699f79df75`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=ShapeShifter` and
`HS_PROFILE_WINDOW=16`; `just profile ShapeShifter` performs the
locked build, flash, capture, and artifact attestation.
