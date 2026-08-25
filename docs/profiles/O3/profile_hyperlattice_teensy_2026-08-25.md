# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-08-25, **global -O3**)

Global-O3 twin of the [shipping report](../shipping/profile_hyperlattice_teensy_2026-08-25.md).
Raw capture: `build/prof/hyperlattice_transition_o3.log`. This replaces the
2026-08-24 two-preset transition report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, two presets, tip `859ab5c03d476` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 170 s full cycle, epoch stretched to 1,600 revolutions |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile_o3 170 16 '-D HS_PROFILE_EPOCH_REVS=1600'` |

Image size: `FLASH: code:70,024, data:147,436, headers:8,844` / `RAM1:
variables:315,040, code:30,920, padding:1,848, free:176,480` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 2,065–2,080 root counter cycles ÷ 600 MHz
matches the measured wall sum within **1.3 ppm**.

## Frame cadence

Peak frame render is **47.79 ms** at frames 2,065–2,080, with **0/2,688
frames spilled**. That leaves 14.71 ms against the 62.5 ms deadline.

The capture covers repeated 320-frame holds and 240-frame transitions in both
directions. Every phase sustains 16 fps. `canvas_buffer_wait` is display-sync
idle, not render work.

## Phase-by-phase readout

### Worst transition window (frames 2,065–2,080)

```text
frame                 62.33 ms  37.40 Mcyc  100%
  hl_shader_draw      41.98 ms  25.19 Mcyc   67%
  hl_timeline_step     9.31 us   5.61 kcyc    0%
  canvas_clear        86.75 us  52.05 kcyc    0%
  canvas_buffer_wait  17.77 ms  10.66 Mcyc   28%
```

Wall min/avg/max = 58.75/62.33/66.86 ms. The exact render peak is 47.79 ms;
the window-average scopes above are attribution only.

### Per-preset table

| # | Preset | Peak render | Spilled | Cadence |
|---:|---|--:|--:|--:|
| 1 | `cubic-flight` | 47.79 ms | 0/1,437 | 16 fps |
| 2 | `hypercube-flight` | 46.17 ms | 0/1,251 | 16 fps |

Initial unlabeled frames are folded into preset 1. The cycle visits both
presets and wraps back to preset 1.

### Per-pixel figures

The draw covers one ≈10,368-sample quadrant. The same two shader bodies are
present as in shipping; global O3 changes their codegen but not dispatch shape.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1149.8/frame  min/avg/max 0.64/1.68/16.71 us  cpu 3.08%
isr_pack         143.7/frame  min/avg/max 5.99/6.49/8.59 us   cpu 1.49%
isr_dma_submit   143.7/frame  min/avg/max 0.77/0.95/9.64 us   cpu 0.21%
```

- Packing remains the dominant DMA-marshaling CPU cost.
- LED transfer is asynchronous; display synchronization is isolated in
  `canvas_buffer_wait`.
- ISR time is included in render scopes because `CYCCNT` free-runs.

## Summary ranking

1. `hl_shader_draw` — dominant render scope; the worst transition retains
   14.71 ms of deadline margin.
2. `canvas_buffer_wait` — display-sync idle and not render work.
3. Column ISR wake and packing — 4.57% aggregate CPU share in the worst window.

Global O3 lowers the shipping peak from 49.67 to 47.79 ms, a **1.88 ms
(3.8%)** improvement. It adds **10,184 bytes** of flash code and **7,936
bytes** of ITCM code, so the selective-O3 shipping image remains the better
full-roster tradeoff.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- No deep per-event scopes were enabled in this timing capture.
- Global O3 is a comparison image, not the shipping configuration.
- The capture was built from clean landed tip `859ab5c03d476`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`; the reproduction command captures the global-O3 twin.
