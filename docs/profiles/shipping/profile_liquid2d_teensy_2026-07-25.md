# Liquid2D on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Liquid2D`).
Raw capture: `build/prof/liquid2d_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The stereographic shader draw (`lq_shader_draw`) is the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Liquid2D 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 70 s capture, 2-preset cycle |
| Reproduce | `bash tools/profile_one.sh Liquid2D profile 70 16` |

Image size: `FLASH: code:47388, data:1062388, headers:8424` / `RAM1:
variables:314336, code:29480, padding:3288, free:177184` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 97–112 root counter cyc ÷ 600 MHz matches
the measured wall sum within **1.2 ppm**.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... buckets`): both presets green — worst
preset peak frame render **29.15 ms** (preset 0), spilled **0/1088 frames (0%)**
— 🟢. Render ~28 ms; `lq_buffer_wait` ~54% — holds 16 fps with ~34 ms of headroom.

## Phase-by-phase readout

Phase schedule: Liquid2D alternates 2 presets (`Preset:` marker). Both are held
regimes with near-identical cost; the worst is below.

### Preset 0 hold (window frames 1073–1088)

```
frame              62.44 ms  37.46 Mcyc  100%
  lq_shader_draw   28.31 ms  16.99 Mcyc   45%
  lq_timeline_step  0.06 ms   0.04 Mcyc    0%
  lq_buffer_wait   34.06 ms  20.44 Mcyc   54%
```

Wall min/avg/max = 61.85/62.44/62.88 ms; percentages are of the parent scope.
`lq_shader_draw` — the per-pixel stereographic shader — is the entire render.

### Per-preset table

Cycle wrap confirmed (`cycle wraps to 0`, both presets visited). Peak/spilled
from `parse_profile.py ... buckets`.

| # | Preset | shader ms/f | peak ms | spilled | fps |
|---|---|--:|--:|--:|--:|
| 0 | preset 0 | 28.31 | 🟢 29.15 | 0/621 (0%) | 16 |
| 1 | preset 1 | ~28 | 🟢 29.06 | 0/467 (0%) | 16 |

### Per-pixel figures

No `filter_blend` in the tree — the stereographic shader writes the canvas
directly. `lq_shader_draw` 16.99 Mcyc ÷ 10,368 = **1,639 cyc per pixel**.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.58/1.73/11.51 us  cpu 3.19%
isr_pack         144/frame  min/avg/max 6.10/6.58/9.46 us  cpu 1.51%
isr_dma_submit   144/frame  min/avg/max 0.72/0.93/1.08 us  cpu 0.21%
```

(window frames 1073–1088.) Total ISR CPU share **4.91%**.

## Summary ranking

1. `lq_shader_draw` — 45% of the frame, 28.3 ms: per-pixel stereographic shader.

Liquid2D holds 16 fps comfortably. The global-`-O3` ceiling buys 29.15 →
26.58 ms (1.10×), no cadence stakes.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- No `filter_blend` — direct per-pixel shade.
- Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Liquid2D`,
`HS_PROFILE_WINDOW=16`; `just profile Liquid2D` = build + flash + capture.
