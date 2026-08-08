# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile SphericalHarmonics`).
Raw capture: `build/prof/sphericalharmonics_ship.log`. Replaces the 2026-07-20
report; half of a matched same-session ship/O3 pair from the full-roster
2026-07-25 sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The field rasterizer (`sh_rasterize`) and terminal blend are the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 220 s capture, epoch stretched to 2048 revs, ordered 24-mode cycle (`HS_PROFILE_ORDERED_CYCLE`). Modes morph back-to-back, so per-mode rows use the peak window nearest each anchor. |
| Reproduce | `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` |

Image size: `FLASH: code:34564, data:535668, headers:8320` / `RAM1:
variables:314368, code:21912, padding:10856, free:177152` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 833–848 root counter cyc ÷ 600 MHz matches
the measured wall sum within **1.6 ppm**.

## Frame cadence

**Pass aggregate** (`buckets`): all 24 modes green — worst peak frame render
**16.13 ms** (mode 20), spilled **0/3488 frames (0%)** — 🟢. Render ~16 ms;
`sh_buffer_wait` ~74% — one of the cheapest cyclers, ~46 ms of headroom.

## Phase-by-phase readout

Phase schedule: SphericalHarmonics morphs continuously through 24 modes
(`Mode:` marker), no held plateau; the regime below is the peak window nearest
the mode-20 anchor.

### Mode 20 peak window (frames 2417–2432)

```
frame            62.26 ms  37.36 Mcyc  100%
  sh_rasterize   15.71 ms   9.43 Mcyc   25%
    filter_blend  1.51 ms   0.91 Mcyc    9%  x10658  85 cyc/blend
  sh_timeline_step 0.02 ms  0.01 Mcyc    0%
  sh_buffer_wait 46.52 ms  27.91 Mcyc   74%
```

Wall min/avg/max = 59.32/62.26/62.53 ms; percentages are of the parent scope.
`sh_rasterize` — the harmonic-field rasterizer — is the whole render, with
`filter_blend` a 9% tail (10,658 blended px = 1.03× the quadrant).

### Per-preset table

Cycle wrap confirmed (`cycle wraps to 0`, all 24 modes visited). Peak/spilled
from `buckets`; peak-window figures per the back-to-back-morph method note. All
24 modes green, peaks span **11.8–16.1 ms**; the four costliest are modes
20 (16.13), 21 (16.11), 22 (15.36), 19 (15.30).

### Per-pixel figures

Quadrant = 10,368 px; `filter_blend` = 10,658 blended px/frame ⇒ **1.03×
coverage**. `filter_blend` 85 cyc/blend; `sh_rasterize` 9.43 Mcyc ÷ 10,658 =
885 cyc per blended pixel.

## Column-ISR / DMA marshaling cost

```
isr_wake        1148/frame  min/avg/max 0.64/1.72/11.51 us  cpu 3.17%
isr_pack         144/frame  min/avg/max 6.11/6.58/9.40 us  cpu 1.51%
isr_dma_submit   144/frame  min/avg/max 0.61/0.93/1.05 us  cpu 0.21%
```

(window frames 2417–2432.) Total ISR CPU share **4.89%**.

## Summary ranking

1. `sh_rasterize` — 25% of the frame, 15.7 ms: harmonic-field rasterizer.
2. `filter_blend` — 1.5 ms terminal blend.

SphericalHarmonics is among the cheapest cyclers; render fits ~4× over in one
window. The global-`-O3` ceiling buys the worst mode 16.13 → 12.85 ms (1.26×).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `sh_rasterize`; its calls ≈ blended px.
- Back-to-back morphs: per-mode rows are peak windows nearest each anchor.
- Dwell knobs (`HS_PROFILE_ORDERED_CYCLE`, `HS_PROFILE_EPOCH_REVS=2048`) change
  how long a mode holds, not per-frame cost. Working tree tip `32478115`, only
  untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=SphericalHarmonics`,
`HS_PROFILE_WINDOW=16`; `just profile SphericalHarmonics` = build + flash +
capture.
