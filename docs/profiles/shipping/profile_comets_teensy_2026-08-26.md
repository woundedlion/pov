# Comets on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Comets`).
Raw capture: `build/prof/comets_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_comets_teensy_2026-07-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses Plot cull/raster and splat `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Comets 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size:

```text
FLASH: code:64,664, data:149,036, headers:8,508
       free for files:1,809,408
RAM1:  variables:315,072, code:35,192, padding:30,344
       free for local variables:143,680
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 3841–3856 root counter
cycles ÷ 600 MHz match the measured wall sum within **3.1 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `cm_draw_trail` averages
16.16 ms/frame; its worst window is 26.52 ms/frame
(frames 3841–3856). Peak frame render is
**33.91 ms** (frames 3841–3856);
spilled **0/4128 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 28.59 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 3841–3856)

```text
frame                         62.23 ms  37.34 Mcyc  100%
  cm_draw_trail               26.52 ms  15.91 Mcyc   43%
    filter_blend               1.38 ms 830.78 kcyc    2%  x13981.9 0.1us/c
  cm_wipe_rebake               2.17 ms   1.30 Mcyc    3%  x1.0 2167.8us/c
  cm_timeline_step             81.2 us  48.78 kcyc    0%  x1.0 81.2us/c
  canvas_clear                 85.3 us  51.19 kcyc    0%  x1.0 85.3us/c
  canvas_buffer_wait          33.37 ms  20.02 Mcyc   54%  x1.0 33372.9us/c
```

Wall min/avg/max = 53.74/62.23/72.17 ms. `cm_draw_trail`
accounts for 42.6% of this measured frame at 26.52 ms/frame.
The complete render is 28.85 ms; `canvas_buffer_wait` contributes
33.37 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–16)

```text
frame                         54.93 ms  32.96 Mcyc  100%
  cm_draw_trail                2.12 ms   1.27 Mcyc    4%
    filter_blend              141.9 us  85.14 kcyc    0%  x1438.4 0.1us/c
  cm_wipe_rebake                0.1 us      49 cyc    0%  x1.0 0.1us/c
  cm_timeline_step             72.9 us  43.77 kcyc    0%  x1.0 72.9us/c
  canvas_clear                 86.4 us  51.87 kcyc    0%  x1.0 86.4us/c
  canvas_buffer_wait          52.65 ms  31.59 Mcyc   96%  x1.0 52646.8us/c
```

Wall min/avg/max = 1.48/54.93/64.09 ms. `cm_draw_trail`
accounts for 3.9% of this measured frame at 2.12 ms/frame.
The complete render is 2.28 ms; `canvas_buffer_wait` contributes
52.65 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `cm_draw_trail` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `cm_draw_trail` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `1` | — | 20/20 | 13,981.9 | 26.52 | 33.91 | 16.0 |
| 2 | `11` | — | 20/20 | 13,647.2 | 24.80 | 28.89 | 16.0 |
| 3 | `10` | — | 20/20 | 12,214.8 | 24.80 | 29.74 | 16.0 |
| 4 | `2` | — | 28/28 | 12,295.1 | 23.96 | 27.75 | 16.0 |
| 5 | `12` | — | 20/20 | 11,761.1 | 23.80 | 29.97 | 16.0 |
| 6 | `4` | — | 20/20 | 12,142.7 | 22.94 | 30.04 | 16.0 |
| 7 | `8` | — | 20/20 | 9,368.2 | 22.31 | 25.25 | 16.0 |
| 8 | `5` | — | 20/20 | 10,613.2 | 22.27 | 29.84 | 16.0 |
| 9 | `9` | — | 20/20 | 9,767.1 | 21.80 | 25.55 | 16.0 |
| 10 | `6` | — | 20/20 | 9,345.7 | 18.80 | 24.55 | 16.0 |
| 11 | `3` | — | 20/20 | 9,510.9 | 18.54 | 22.97 | 16.0 |
| 12 | `0` | — | 10/10 | 8,353.9 | 14.15 | 16.67 | 16.0 |
| 13 | `7` | — | 20/20 | 5,837.6 | 12.95 | 17.13 | 16.0 |

### Per-pixel figures

The peak dominant-scope window blends 13,981.9 pixels/frame, 1.35× the 10,368-pixel quadrant. `filter_blend` costs 59.4 cycles/blend; `cm_draw_trail` uses 1138.0 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1153.8/f  0.52/1.59/19.35 us  cpu 2.94%
isr_pack          143.9/f  6.23/6.74/13.44 us  cpu 1.55%
isr_dma_submit    143.9/f  0.58/0.94/11.95 us  cpu 0.22%
```

- Pack plus submit costs 1.10 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.71% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `cm_draw_trail` — 25.9% of aggregated root time, 16.16 ms/frame: inclusive measured scope in the live driver.
2. `filter_blend` — 1.1% of aggregated root time, 0.69 ms/frame: inclusive measured scope in the live driver.
3. `cm_wipe_rebake` — 1.0% of aggregated root time, 0.63 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches Plot cull/raster and splat `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=2400`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=Comets` and
`HS_PROFILE_WINDOW=16`; `just profile Comets` performs the
locked build, flash, capture, and artifact attestation.
