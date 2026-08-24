# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-08-24, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Raw capture: `build/prof/hyperlattice_ship.log`. First profile of this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, four presets, tip `374d1a2d6391` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 380 s capture; epoch stretched to 3,600 revolutions |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile 380 16 '-D HS_PROFILE_EPOCH_REVS=3600'` |

Image size: `FLASH: code:53,200, data:147,132, headers:8,564` / `RAM1:
variables:315,008, code:21,976, padding:10,792, free:176,512` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 1,649–1,664 root counter cycles ÷ 600 MHz
match the measured wall sum within **1.5 ppm**.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`hl_shader_draw` reaches 167.43 ms/f at the slowest clean preset window; peak
frame render is **172.06 ms** at frames 1,665–1,680, with **2,672/2,672**
frames spilled (100%).

A display window is 62.5 ms. HyperLattice renders one quadrant ≈ 10,368
samples per frame. Every captured preset exceeds the display budget, with
worst-window cadence falling to 5.3 fps.

## Phase-by-phase readout

The 3,600-revolution profile epoch permits the default 320-frame dwell and
240-frame segue to complete a full four-preset wrap. The parser's unlabeled
startup bucket was folded into preset 1; the preset marker then confirms the
same `cubic-flight` phase before advancing through presets 2–4 and wrapping.

### Worst window (frames 1,665–1,680)

```text
frame                 186.35 ms 111.81 Mcyc  100%
  hl_shader_draw      163.43 ms  98.06 Mcyc   87%
  fx_timeline_step     12.4 us    7.46 kcyc    0%
  canvas_clear         85.5 us   51.33 kcyc    0%
  canvas_buffer_wait   20.34 ms  12.20 Mcyc   10%
```

Wall min/avg/max = 175.97/186.35/192.51 ms. The layered reflected-lattice
shader owns 87% of wall time. Display-sync wait reflects the accumulated
missed cadence rather than usable render margin.

### Preset sweep

| # | Preset | Clean shader max | Peak render | Spilled | Worst cadence |
|---:|---|---:|---:|---:|---:|
| 1 | `cubic-flight` | 161.49 ms | 171.72 ms | 878/878 | 5.3 fps |
| 2 | `deep-grid` | 119.06 ms | 122.13 ms | 676/676 | 8 fps |
| 3 | `dimensional-rift` | 139.83 ms | 152.80 ms | 559/559 | 5.3 fps |
| 4 | `hypercube-flight` | 167.43 ms | 172.06 ms | 559/559 | 5.3 fps |

### Per-pixel figures

The draw covers one 10,368-sample quadrant. The slowest clean preset window
spends about 9,690 cycles per sample in `hl_shader_draw`; no separate blend
scope is enabled, so this includes the whole layered shader evaluation.

## Column-ISR / DMA marshaling cost

```text
isr_wake       3435.9/frame  min/avg/max 0.80/2.02/14.87 us  cpu 3.72%
isr_pack        429.4/frame  min/avg/max 6.25/6.93/9.43 us   cpu 1.59%
isr_dma_submit  429.4/frame  min/avg/max 0.77/0.94/6.85 us   cpu 0.21%
```

- Packing dominates submit CPU cost; wire transfer remains asynchronous.
- The measured ISR share is included in the free-running cycle scopes.
- The shader alone exceeds the 62.5 ms cadence budget in every preset.

## Summary ranking

1. `hl_shader_draw` — 87% of the worst window, 163.43 ms/frame.
2. `canvas_buffer_wait` — 10%, 20.34 ms/frame.
3. Column ISR wake and packing — 5.31% aggregate CPU share.

The global-O3 twin peaks at 148.24 ms, 13.8% below the shipping peak, but
still spills every frame.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- No per-pixel or blend scope was enabled, avoiding instrumentation distortion.
- The shipping image uses the repository's selective `HS_O3` regions on an
  otherwise `-Os` build.
- The extended epoch changes only profiling cadence, not the effect preset
  definitions or render path.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=3600`; the reproduction command
above builds, flashes and captures the full preset wrap.
