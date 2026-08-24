# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-08-24, **global -O3**)

Point-in-time snapshot (regenerate with `just profile-o3 HyperLattice`).
Raw capture: `build/prof/hyperlattice_o3.log`. Global-O3 twin of the shipping
capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3` env: global `-O3` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, four presets, tip `374d1a2d6391` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 380 s capture; epoch stretched to 3,600 revolutions |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile_o3 380 16 '-D HS_PROFILE_EPOCH_REVS=3600'` |

Image size: `FLASH: code:64,112, data:147,248, headers:8,800` / `RAM1:
variables:315,008, code:30,616, padding:2,152, free:176,512` / `RAM2:
variables:520,064, free:4,224`. Relative to shipping, global `-O3` adds
10,912 bytes of flash code and 8,640 bytes of ITCM code.

Exactness cross-check: window frames 1,649–1,664 root counter cycles ÷ 600 MHz
match the measured wall sum within **0.7 ppm**.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`hl_shader_draw` reaches 143.80 ms/f at the slowest clean preset window; peak
frame render is **148.24 ms** at frames 1,649–1,664, with **2,816/2,816**
frames spilled (100%).

A display window is 62.5 ms. Global `-O3` reduces every preset's shader cost,
but no captured preset reaches the display budget. Worst-window cadence remains
5.3 fps.

## Phase-by-phase readout

The 3,600-revolution profile epoch permits the default 320-frame dwell and
240-frame segue to complete a full four-preset wrap. The parser's unlabeled
startup bucket was folded into preset 1, matching the shipping report.

### Worst window (frames 1,649–1,664)

```text
frame                 187.50 ms 112.50 Mcyc  100%
  hl_shader_draw      143.80 ms  86.28 Mcyc   76%
  fx_timeline_step     10.5 us    6.28 kcyc    0%
  canvas_clear         84.7 us   50.84 kcyc    0%
  canvas_buffer_wait   40.89 ms  24.53 Mcyc   21%
```

Wall min/avg/max = 184.86/187.50/190.50 ms. The layered reflected-lattice
shader owns 76% of wall time. Display-sync wait reflects accumulated missed
cadence rather than usable render margin.

### Preset sweep

| # | Preset | Clean shader max | Peak render | Spilled | Worst cadence |
|---:|---|---:|---:|---:|---:|
| 1 | `cubic-flight` | 138.46 ms | 148.11 ms | 878/878 | 5.3 fps |
| 2 | `deep-grid` | 97.23 ms | 100.88 ms | 820/820 | 8 fps |
| 3 | `dimensional-rift` | 117.65 ms | 126.21 ms | 559/559 | 5.3 fps |
| 4 | `hypercube-flight` | 143.80 ms | 148.24 ms | 559/559 | 5.3 fps |

### Per-pixel figures

The draw covers one 10,368-sample quadrant. The slowest clean preset window
spends about 8,322 cycles per sample in `hl_shader_draw`; no separate blend
scope is enabled, so this includes the whole layered shader evaluation.

## Column-ISR / DMA marshaling cost

```text
isr_wake       3456.9/frame  min/avg/max 0.61/1.84/15.55 us  cpu 3.38%
isr_pack        432.1/frame  min/avg/max 6.00/6.77/9.05 us   cpu 1.55%
isr_dma_submit  432.1/frame  min/avg/max 0.78/0.94/7.76 us   cpu 0.21%
```

- Packing dominates submit CPU cost; wire transfer remains asynchronous.
- The measured ISR share is included in the free-running cycle scopes.
- The shader alone exceeds the 62.5 ms cadence budget in every preset.

## Summary ranking

1. `hl_shader_draw` — 76% of the worst window, 143.80 ms/frame.
2. `canvas_buffer_wait` — 21%, 40.89 ms/frame.
3. Column ISR wake and packing — 4.93% aggregate CPU share.

Global `-O3` lowers the peak from 172.06 ms to 148.24 ms (13.8%), at a cost
of 10,912 bytes of flash code and 8,640 bytes of ITCM code. It does not change
the all-frame spill result.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- No per-pixel or blend scope was enabled, avoiding instrumentation distortion.
- This is a global `-O3` comparison build, not the shipping configuration.
- The extended epoch changes only profiling cadence, not the effect preset
  definitions or render path.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS=3600`; the reproduction command
above builds, flashes and captures the full preset wrap.
