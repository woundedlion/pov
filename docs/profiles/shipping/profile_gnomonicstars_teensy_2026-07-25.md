# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GnomonicStars`).
Raw capture: `build/prof/gnomonicstars_ship.log`. Replaces the 2026-07-20
report; half of a matched same-session ship/O3 pair from the full-roster
2026-07-25 sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The star scan (`gn_star_scan`) and terminal blend are the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GnomonicStars 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh GnomonicStars profile 70 32` |

Image size: `FLASH: code:44052, data:537780, headers:9008` / `RAM1:
variables:314432, code:29384, padding:3384, free:177088` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 353–384 root counter cyc ÷ 600 MHz matches
the measured wall sum within **2.6 ppm**.

## Frame cadence

**Pass aggregate**: peak frame render **38.07 ms** (frames 353–384), spilled
**0/1088 frames (0%)** — 🟢. Render averages ~24 ms with peaks near 38 ms;
`gn_buffer_wait` is 61% of the frame, so it holds 16 fps with ~24 ms of headroom.

## Phase-by-phase readout

Phase schedule: GnomonicStars continuously rasterizes a rotating star field;
single regime.

### Peak window (frames 353–384)

```
frame              62.92 ms  37.75 Mcyc  100%
  gn_draw_stars    24.29 ms  14.58 Mcyc   38%
    gn_star_scan   23.41 ms  14.05 Mcyc   96%  x600.0  39 us/star
      filter_blend  1.12 ms   0.67 Mcyc    4%  x6748  99 cyc/blend
  gn_timeline_step  0.06 ms   0.03 Mcyc    0%
  gn_buffer_wait   38.58 ms  23.15 Mcyc   61%
```

Wall min/avg/max = 48.03/62.92/78.28 ms; percentages are of the parent scope.
`gn_star_scan` is 96% of the star draw — 600 star scans at 39 us each — with a
thin `filter_blend` tail (6,748 blended px = 0.65× the quadrant, so coverage is
sparse).

### Per-pixel figures

Quadrant = 10,368 px; `filter_blend` = 6,748 blended px/frame ⇒ **0.65×
coverage** (sparse star field). `filter_blend` 99 cyc/blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1160/frame  min/avg/max 0.47/1.74/13.78 us  cpu 3.21%
isr_pack         145/frame  min/avg/max 6.11/6.56/9.41 us  cpu 1.51%
isr_dma_submit   145/frame  min/avg/max 0.60/0.93/4.15 us  cpu 0.21%
```

(window frames 353–384.) Total ISR CPU share **4.93%**.

## Summary ranking

1. `gn_star_scan` — 96% of the star draw, 23.4 ms: 600 star scans.
2. `filter_blend` — 1.1 ms terminal blend.

GnomonicStars holds 16 fps comfortably (~24 ms margin). The global-`-O3` ceiling
buys 38.07 → 32.47 ms (1.17×) — a healthy codegen win with no cadence stakes.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `gn_star_scan`; its calls ≈ blended px.
- Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=GnomonicStars`,
`HS_PROFILE_WINDOW=32`; `just profile GnomonicStars` = build + flash + capture.
