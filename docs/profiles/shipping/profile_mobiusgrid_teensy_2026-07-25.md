# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MobiusGrid`).
Raw capture: `build/prof/mobiusgrid_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The curve raster (`mg_draw_grid`) is the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MobiusGrid 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh MobiusGrid profile 70 32` |

Image size: `FLASH: code:60332, data:1062188, headers:8992` / `RAM1:
variables:314432, code:41384, padding:24152, free:144320` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 225–256 root counter cyc ÷ 600 MHz matches
the measured wall sum within **3.6 ppm**.

## Frame cadence

**Pass aggregate**: peak frame render **17.42 ms** (frames 609–640), spilled
**0/1088 frames (0%)** — 🟢. Render ~13 ms; `mg_buffer_wait` 78% — holds 16 fps
with ~45 ms of headroom.

## Phase-by-phase readout

Phase schedule: MobiusGrid continuously draws a warped line/ring grid; single
regime.

### Peak window (frames 609–640)

```
frame              62.50 ms  37.50 Mcyc  100%
  mg_draw_grid     12.77 ms   7.66 Mcyc   20%
    mg_lines_draw   7.75 ms   4.65 Mcyc   60%
    mg_rings_draw   5.02 ms   3.01 Mcyc   39%
  mg_wipe_rebake    0.38 ms   0.23 Mcyc    0%
  mg_timeline_step  0.15 ms   0.09 Mcyc    0%
  mg_buffer_wait   49.19 ms  29.51 Mcyc   78%
```

Wall min/avg/max = 54.25/62.50/70.66 ms; percentages are of the parent scope.
The curve raster splits ~60/40 between line and ring spans; both are cheap.

### Per-pixel figures

No `filter_blend` in the tree — the curve rasterizer writes spans directly.
Render is curve-length-bound.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.64/1.72/15.94 us  cpu 3.16%
isr_pack         144/frame  min/avg/max 6.11/6.57/9.47 us  cpu 1.51%
isr_dma_submit   144/frame  min/avg/max 0.70/0.92/5.83 us  cpu 0.21%
```

(window frames 609–640.) Total ISR CPU share **4.88%**.

## Summary ranking

1. `mg_lines_draw` — 7.8 ms: warped line spans.
2. `mg_rings_draw` — 5.0 ms: ring spans.

MobiusGrid is among the cheapest effects; render fits ~3.5× over in one window.
The global-`-O3` ceiling buys 17.42 → 15.36 ms (1.13×), immaterial to cadence.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs). No `filter_blend` — direct span
  raster. Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MobiusGrid`,
`HS_PROFILE_WINDOW=32`; `just profile MobiusGrid` = build + flash + capture.
