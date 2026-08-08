# RingShower on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingShower`).
Raw capture: `build/prof/ringshower_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The ring plot (`rsh_ring_plot`) is the hot path. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingShower 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingShower profile 70 32` |

Image size: `FLASH: code:52668, data:1060756, headers:8872` / `RAM1:
variables:314368, code:36936, padding:28600, free:144384` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 961–992 root counter cyc ÷ 600 MHz matches
the measured wall sum within **0.1 ppm**.

## Frame cadence

**Pass aggregate**: peak frame render **4.17 ms** (frames 993–1024), spilled
**0/1088 frames (0%)** — 🟢. **The cheapest effect on the roster**:
`rsh_buffer_wait` is 96% of the frame, so render fits ~15× over in one display
window.

## Phase-by-phase readout

Phase schedule: RingShower drops a sparse stream of rings; single regime.

### Peak window (frames 993–1024)

```
frame               62.43 ms  37.46 Mcyc  100%
  rsh_draw_rings     2.44 ms   1.46 Mcyc    3%
    rsh_ring_plot    2.44 ms   1.46 Mcyc   99%  x3.3  740 us/ring
  rsh_timeline_step  0.02 ms   0.01 Mcyc    0%
  rsh_buffer_wait   59.97 ms  35.98 Mcyc   96%
```

Wall min/avg/max = 60.68/62.43/63.75 ms; percentages are of the parent scope.
Only ~3.3 rings are plotted per frame; the effect is nearly all idle.

### Per-pixel figures

No `filter_blend` in the tree — the ring plot writes spans directly. Render is
ring-count-bound (~3.3 rings/frame).

## Column-ISR / DMA marshaling cost

```
isr_wake        1151/frame  min/avg/max 0.64/1.70/11.44 us  cpu 3.12%
isr_pack         144/frame  min/avg/max 6.11/6.55/9.50 us  cpu 1.50%
isr_dma_submit   144/frame  min/avg/max 0.60/0.93/1.07 us  cpu 0.21%
```

(window frames 993–1024.) Total ISR CPU share **4.83%** — larger than the 4.17 ms
render itself, since the column marshaling is a fixed per-frame cost regardless
of how little the effect draws.

## Summary ranking

1. `rsh_ring_plot` — 2.4 ms: ~3.3 ring plots. The entire (tiny) render.

RingShower is the roster's cheapest effect by a wide margin. The global-`-O3`
ceiling buys 4.17 → 3.93 ms — irrelevant.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs); here the ISR cost exceeds the
  render cost.
- No `filter_blend` — direct span raster. Working tree tip `32478115`, only
  untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingShower`,
`HS_PROFILE_WINDOW=32`; `just profile RingShower` = build + flash + capture.
