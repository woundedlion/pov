# Comets on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Comets`).
Raw capture: `build/prof/comets_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The point/trail scan (`cm_point_scan`) is the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Comets 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, epoch stretched to 2400 revs, 12-preset cycle |
| Reproduce | `bash tools/profile_one.sh Comets profile 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size: `FLASH: code:54428, data:1063240, headers:8724` / `RAM1:
variables:314432, code:35752, padding:29784, free:144320` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 3985–4000 root counter cyc ÷ 600 MHz matches
the measured wall sum within **3.8 ppm**.

## Frame cadence

**Pass aggregate** (`buckets`): all 12 presets green — worst preset peak frame
render **36.27 ms** (preset 0, transition-inclusive), spilled **0/4128 frames
(0%)** — 🟢. Clean hold ~20 ms; `cm_buffer_wait` ~75%.

## Phase-by-phase readout

Phase schedule: Comets steps through 12 presets. The costliest, preset 0, is
below at its clean hold; the bucket peak of 36.27 ms includes its transition.

### Preset 0 clean hold (window frames 3985–4000)

```
frame              63.13 ms  37.88 Mcyc  100%
  cm_draw_trail    15.13 ms   9.08 Mcyc   23%
    cm_point_scan  14.93 ms   8.96 Mcyc   98%  x560.8  27 us/point
      filter_blend  0.98 ms   0.59 Mcyc    6%  x6694  88 cyc/blend
  cm_wipe_rebake    0.00 ms   0.00 Mcyc    0%
  cm_timeline_step  0.09 ms   0.05 Mcyc    0%
  cm_buffer_wait   47.91 ms  28.75 Mcyc   75%
```

Wall min/avg/max = 53.76/63.13/72.78 ms; percentages are of the parent scope.
`cm_point_scan` is 98% of the trail draw — 560.8 points at 27 us each — with a
thin `filter_blend` tail (6,694 blended px = 0.65× the quadrant).

### Per-preset table

Cycle wrap confirmed (`cycle wraps to 0`, all 12 presets visited). Peak/spilled
from `buckets`.

| # | Preset | peak ms | spilled | | # | Preset | peak ms | spilled |
|---|---|--:|--:|---|---|---|--:|--:|
| 0 | 🟢 36.27 | 0/480 | | 4 | 🟢 32.74 | 0/320 |
| 10 | 🟢 35.94 | 0/320 | | 3 | 🟢 32.46 | 0/320 |
| 9 | 🟢 33.30 | 0/320 | | 8 | 🟢 31.17 | 0/320 |
| 7 | 🟢 33.05 | 0/320 | | 1 | 🟢 26.77 | 0/448 |
| 11 | 🟢 32.87 | 0/320 | | 2 | 🟢 26.64 | 0/320 |
| | | | | 5 | 🟢 25.03 | 0/320 |
| | | | | 6 | 🟢 23.18 | 0/320 |

All 12 green; peaks span 23.2–36.3 ms.

### Per-pixel figures

Quadrant = 10,368 px; `filter_blend` = 6,694 blended px/frame ⇒ **0.65×
coverage**. `cm_point_scan` 8.96 Mcyc ÷ 560.8 points = **16 kcyc per point**.

## Column-ISR / DMA marshaling cost

```
isr_wake        1164/frame  min/avg/max 0.64/1.72/11.40 us  cpu 3.16%
isr_pack         146/frame  min/avg/max 6.11/6.55/9.44 us  cpu 1.51%
isr_dma_submit   146/frame  min/avg/max 0.62/0.93/1.02 us  cpu 0.21%
```

(window frames 3985–4000.) Total ISR CPU share **4.88%**.

## Summary ranking

1. `cm_point_scan` — 98% of the trail draw, 14.9 ms: 560.8 point scans.

Comets is a cheap, comfortably-green effect (all 12 presets ≤ 36 ms). The
global-`-O3` ceiling buys the worst preset 36.27 → 30.46 ms.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `cm_point_scan`; its calls ≈ blended px.
- Dwell knob `-D HS_PROFILE_EPOCH_REVS=2400` stretches the epoch, not per-frame
  cost. Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Comets`,
`HS_PROFILE_WINDOW=16`; `just profile Comets` = build + flash + capture.
