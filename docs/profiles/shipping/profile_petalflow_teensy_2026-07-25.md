# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile PetalFlow`).
Raw capture: `build/prof/petalflow_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

The roster index ([`../README.md`](../README.md)) ranks this row from a
2026-07-26 11:37 log that postdates this report: the index peak is the current
figure and the numbers below are the earlier capture. `just profile PetalFlow`
regenerates this report against a fresh log.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The ring scan (`pf_ring_scan`) is the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh PetalFlow profile 70 32` |

Image size: `FLASH: code:48596, data:536136, headers:9180` / `RAM1:
variables:314400, code:34744, padding:30792, free:144352` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 289–320 root counter cyc ÷ 600 MHz matches
the measured wall sum within **2.4 ppm**.

## Frame cadence

**Pass aggregate**: peak frame render **11.81 ms** (frames 289–320), spilled
**0/1088 frames (0%)** — 🟢. Render ~11 ms; `pf_buffer_wait` 82% — the second
cheapest effect on the roster, ~51 ms of headroom.

## Phase-by-phase readout

Phase schedule: PetalFlow continuously draws a small stack of flowing rings;
single regime.

### Peak window (frames 289–320)

```
frame              62.50 ms  37.50 Mcyc  100%
  pf_draw_rings    10.96 ms   6.58 Mcyc   17%
    pf_ring_scan    9.87 ms   5.92 Mcyc   90%  x22.7  435 us/ring
    pf_ring_build   0.98 ms   0.59 Mcyc    8%  x22.7
  pf_timeline_step  0.01 ms   0.01 Mcyc    0%
  pf_buffer_wait   51.52 ms  30.91 Mcyc   82%
```

Wall min/avg/max = 61.59/62.50/63.39 ms; percentages are of the parent scope.
`pf_ring_scan` is 90% of the ring draw — 22.7 rings at 435 us each.

### Per-pixel figures

No `filter_blend` in the tree — the ring scan writes spans directly. Render is
ring-count-bound.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.54/1.70/15.75 us  cpu 3.13%
isr_pack         144/frame  min/avg/max 6.11/6.56/9.48 us  cpu 1.51%
isr_dma_submit   144/frame  min/avg/max 0.64/0.93/8.73 us  cpu 0.21%
```

(window frames 289–320.) Total ISR CPU share **4.85%**.

## Summary ranking

1. `pf_ring_scan` — 90% of the ring draw, 9.9 ms: 22.7 ring scans.
2. `pf_ring_build` — 1.0 ms.

PetalFlow is the second cheapest roster effect; render fits ~5× over in one
window. The global-`-O3` ceiling buys 11.81 → 10.85 ms, immaterial.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs). No `filter_blend` — direct span
  raster. Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=PetalFlow`,
`HS_PROFILE_WINDOW=32`; `just profile PetalFlow` = build + flash + capture.
