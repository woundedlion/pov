# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Voronoi`).
Raw capture: `build/prof/voronoi_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The block-union top-2 shade (`vo_shade`) is the hot path; per the Voronoi resolution it holds 16 fps with no HS_O3 region of its own. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Voronoi 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Voronoi profile 70 32` |

Image size: `FLASH: code:34244, data:535744, headers:8564` / `RAM1:
variables:314368, code:20952, padding:11816, free:177152` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 1025–1056 root counter cyc ÷ 600 MHz matches
the measured wall sum within **3.9 ppm**.

## Frame cadence

**Pass aggregate**: peak frame render **10.27 ms** (frames 577–608), spilled
**0/1088 frames (0%)** — 🟢. Render ~10 ms; `vo_buffer_wait` 83% — holds 16 fps
with ~52 ms of headroom.

## Phase-by-phase readout

Phase schedule: Voronoi continuously shades a moving cell diagram; single regime.

### Peak window (frames 577–608)

```
frame             62.47 ms  37.48 Mcyc  100%
  vo_shade         9.55 ms   5.73 Mcyc   15%
  vo_kdtree        0.29 ms   0.17 Mcyc    0%
  vo_animate       0.16 ms   0.10 Mcyc    0%
  vo_buffer_wait  52.47 ms  31.48 Mcyc   83%
```

Wall min/avg/max = 61.95/62.47/62.90 ms; percentages are of the parent scope.
`vo_shade` — the block-union top-2 nearest-site shade — is the whole render;
`vo_kdtree` (site acceleration) and `vo_animate` are rounding costs.

### Per-pixel figures

No `filter_blend` in the tree — the shade writes the canvas directly.
`vo_shade` 5.73 Mcyc ÷ 10,368 = **553 cyc per pixel** (block-union top-2 site
test per pixel).

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.47/1.73/11.25 us  cpu 3.19%
isr_pack         144/frame  min/avg/max 6.11/6.58/9.37 us  cpu 1.51%
isr_dma_submit   144/frame  min/avg/max 0.63/0.93/1.13 us  cpu 0.21%
```

(window frames 577–608.) Total ISR CPU share **4.91%**.

## Summary ranking

1. `vo_shade` — 15% of the frame, 9.6 ms: block-union top-2 per-pixel shade.
2. `vo_kdtree` / `vo_animate` — 0.29 / 0.16 ms, negligible.

Voronoi is among the cheapest effects (~52 ms margin). The block-union shade was
the landing that took it to a hard 16 fps lock; the global-`-O3` ceiling buys
10.27 → 7.57 ms (1.36×), a large relative win with no cadence stakes.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- No `filter_blend` — direct per-pixel shade.
- Selective -O3: no HS_O3 region touches this effect (per the Voronoi block-union
  resolution); the ship/O3 gap is pure global `-O3` + `-ffast-math`.
- Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Voronoi`,
`HS_PROFILE_WINDOW=32`; `just profile Voronoi` = build + flash + capture.
