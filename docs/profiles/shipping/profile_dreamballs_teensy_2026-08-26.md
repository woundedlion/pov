# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DreamBalls`).
Raw capture: `build/prof/dreamballs_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_dreamballs_teensy_2026-08-25.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses mesh/shape scan, transformer, and raster `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DreamBalls 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 230 s capture, `-D HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` |

Image size:

```text
FLASH: code:110,304, data:186,944, headers:8,928
       free for files:1,725,440
RAM1:  variables:315,328, code:50,008, padding:15,528
       free for local variables:143,424
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2849–2864 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.8 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `db_timeline_step` averages
23.48 ms/frame; its worst window is 37.29 ms/frame
(frames 2849–2864). Peak frame render is
**38.73 ms** (frames 2593–2608);
spilled **0/3648 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 23.77 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2849–2864)

```text
frame                         62.62 ms  37.57 Mcyc  100%
  db_timeline_step            37.29 ms  22.38 Mcyc   60%
    db_draw                   37.24 ms  22.35 Mcyc   59%
      db_draw_scene           37.24 ms  22.35 Mcyc   59%
        db_mesh_plot          35.80 ms  21.48 Mcyc   57%
          filter_blend         5.54 ms   3.32 Mcyc    9%  x41936.6 0.1us/c
        db_orient             255.7 us 153.44 kcyc    0%  x6.0 42.6us/c
        db_displace            1.03 ms 615.90 kcyc    2%  x6.0 171.1us/c
  canvas_clear                 86.6 us  52.00 kcyc    0%  x1.0 86.6us/c
  canvas_buffer_wait          25.24 ms  15.14 Mcyc   40%  x1.0 25236.7us/c
```

Wall min/avg/max = 59.26/62.62/65.87 ms. `db_timeline_step`
accounts for 59.6% of this measured frame at 37.29 ms/frame.
The complete render is 37.38 ms; `canvas_buffer_wait` contributes
25.24 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1585–1600)

```text
frame                         61.68 ms  37.01 Mcyc  100%
  db_timeline_step            11.55 ms   6.93 Mcyc   19%
    db_draw                   11.48 ms   6.89 Mcyc   19%
      db_draw_scene           11.48 ms   6.89 Mcyc   19%
        db_mesh_plot          11.15 ms   6.69 Mcyc   18%
          filter_blend         1.81 ms   1.08 Mcyc    3%  x13660.8 0.1us/c
        db_orient              53.6 us  32.17 kcyc    0%  x3.8 14.3us/c
        db_displace           211.3 us 126.82 kcyc    0%  x3.8 56.4us/c
  canvas_clear                 84.7 us  50.82 kcyc    0%  x1.0 84.7us/c
  canvas_buffer_wait          50.04 ms  30.02 Mcyc   81%  x1.0 50040.7us/c
```

Wall min/avg/max = 50.55/61.68/63.49 ms. `db_timeline_step`
accounts for 18.7% of this measured frame at 11.55 ms/frame.
The complete render is 11.64 ms; `canvas_buffer_wait` contributes
50.04 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `db_timeline_step` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `db_timeline_step` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `9` | — | 20/20 | 41,936.6 | 37.29 | 38.73 | 16.0 |
| 2 | `1` | — | 20/20 | 42,649.6 | 35.59 | 36.55 | 16.0 |
| 3 | `0` | — | 20/20 | 42,451.1 | 35.19 | 36.07 | 16.0 |
| 4 | `8` | — | 20/20 | 32,549.9 | 27.76 | 29.36 | 16.0 |
| 5 | `10` | — | 20/20 | 32,737.9 | 27.30 | 29.49 | 16.0 |
| 6 | `4` | — | 20/20 | 26,032.2 | 21.54 | 24.18 | 16.0 |
| 7 | `3` | — | 20/20 | 23,171.1 | 19.97 | 21.64 | 16.0 |
| 8 | `2` | — | 28/28 | 23,130.0 | 19.62 | 20.48 | 16.0 |
| 9 | `7` | — | 20/20 | 21,331.6 | 18.37 | 19.18 | 16.0 |
| 10 | `6` | — | 20/20 | 14,586.1 | 12.88 | 14.16 | 16.0 |
| 11 | `5` | — | 20/20 | 14,874.1 | 12.62 | 13.21 | 16.0 |

### Per-pixel figures

The peak dominant-scope window blends 41,936.6 pixels/frame, 4.04× the 10,368-pixel quadrant. `filter_blend` costs 79.2 cycles/blend; `db_timeline_step` uses 533.6 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.1/f  0.51/1.68/13.52 us  cpu 3.10%
isr_pack          143.9/f  6.23/6.83/9.65 us  cpu 1.57%
isr_dma_submit    143.9/f  0.60/0.95/2.64 us  cpu 0.22%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.89% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `db_timeline_step` — 37.6% of aggregated root time, 23.48 ms/frame: inclusive measured scope in the live driver.
2. `db_draw` — 37.5% of aggregated root time, 23.42 ms/frame: inclusive measured scope in the live driver.
3. `db_draw_scene` — 37.5% of aggregated root time, 23.42 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches mesh/shape scan, transformer, and raster `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=2000`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=DreamBalls` and
`HS_PROFILE_WINDOW=16`; `just profile DreamBalls` performs the
locked build, flash, capture, and artifact attestation.
