# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-08-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`). Raw
capture: `build/prof/hyperlattice_transition_ship.log`. This replaces the
2026-08-24 two-preset transition report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with `HS_O3` on the cached shader and specialized slice hot path |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, two presets, tip `859ab5c03d476` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 170 s full cycle, epoch stretched to 1,600 revolutions |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile 170 16 '-D HS_PROFILE_EPOCH_REVS=1600'` |

Image size: `FLASH: code:59,840, data:147,248, headers:8,976` / `RAM1:
variables:315,008, code:22,984, padding:9,784, free:176,512` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 2,065–2,080 root counter cycles ÷ 600 MHz
matches the measured wall sum within **4.5 ppm**.

## Frame cadence

Peak frame render is **49.67 ms** at frames 2,065–2,080, with **0/2,688
frames spilled**. That leaves 12.83 ms against the 62.5 ms display deadline
and clears the original 59 ms goal by 9.33 ms.

The capture covers repeated 320-frame holds and 240-frame transitions in both
directions. Every hold and every transition frame sustains 16 fps. The
`canvas_buffer_wait` scope is display-sync idle, not render work.

## Phase-by-phase readout

### Worst transition window (frames 2,065–2,080)

```text
frame                 62.33 ms  37.40 Mcyc  100%
  hl_shader_draw      43.68 ms  26.21 Mcyc   70%
  hl_timeline_step    15.38 us   9.26 kcyc    0%
  canvas_clear        87.63 us  52.59 kcyc    0%
  canvas_buffer_wait  16.07 ms   9.64 Mcyc   25%
```

Wall min/avg/max = 58.22/62.33/67.20 ms. The exact render peak is 49.67 ms;
the window-average scopes above are retained only for hot-path attribution.

### Per-preset table

| # | Preset | Peak render | Spilled | Cadence |
|---:|---|--:|--:|--:|
| 1 | `cubic-flight` | 49.67 ms | 0/1,437 | 16 fps |
| 2 | `hypercube-flight` | 47.71 ms | 0/1,251 | 16 fps |

Initial unlabeled frames are folded into preset 1. The cycle visits both
presets and wraps back to preset 1.

### Per-pixel figures

The draw covers one ≈10,368-sample quadrant. ARM codegen contains only the
existing generic flash shader and one specialized ITCM shader; the transition
dispatch does not instantiate a third body.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1149.9/frame  min/avg/max 0.67/1.80/16.29 us  cpu 3.31%
isr_pack         143.8/frame  min/avg/max 6.24/6.69/9.00 us   cpu 1.54%
isr_dma_submit   143.8/frame  min/avg/max 0.81/0.95/8.82 us   cpu 0.21%
```

- Packing remains the dominant DMA-marshaling CPU cost.
- LED transfer is asynchronous; display synchronization is isolated in
  `canvas_buffer_wait`.
- ISR time is included in render scopes because `CYCCNT` free-runs.

## Summary ranking

1. `hl_shader_draw` — dominant render scope; its worst transition still leaves
   12.83 ms of measured deadline margin.
2. `canvas_buffer_wait` — display-sync idle and not optimization headroom.
3. Column ISR wake and packing — 4.85% aggregate CPU share in the worst window.

The previous transition-inclusive shipping capture peaked at 65.37 ms and
spilled 104/2,592 frames. Dispatching all Chrome/Depth/two-shell 4D frames
through the ITCM specialization reduces peak by **15.70 ms (24.0%)** and
eliminates every spill.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- No deep per-event scopes were enabled in this timing capture.
- The specialized distance/fog-AA path can differ only in faint fringe
  coverage; regression tests bound premultiplied visible error.
- The capture and installed WASM artifact were built from clean landed tip
  `859ab5c03d476`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`; the reproduction command builds, flashes, and captures
the complete two-preset cycle.
