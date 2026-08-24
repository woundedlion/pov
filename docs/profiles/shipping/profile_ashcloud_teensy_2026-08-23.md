# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-08-23, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `build/prof/ashcloud_ship.log`. First profile of this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions active |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, tip `b203edd365fc` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AshCloud profile 70 32` |

Image size: `FLASH: code:69,896, data:151,504, headers:9,000` / `RAM1:
variables:315,072, code:26,120, padding:6,648, free:176,448` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 545–576 root counter cycles ÷ 600 MHz
match the measured wall sum within **2.7 ppm**.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`fx_shader_draw` reaches 45.45 ms/f at the slowest window; peak frame render
is **52.39 ms** at frames 513–544, with **0/1,088** frames spilled (0%).

A display window is 62.5 ms. AshCloud renders one quadrant ≈ 10,368 samples
per frame and keeps 10.11 ms of worst-frame margin, so it holds 16 fps for
the whole capture. `canvas_buffer_wait` is the round-up idle to the next
display flip.

## Phase-by-phase readout

The single preset evolves its noise field continuously; frames 513–544 are
the capture's worst render window.

### Worst window (frames 513–544)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  fx_shader_draw       45.20 ms  27.12 Mcyc   72%
  fx_prepare_frame      3.87 ms   2.32 Mcyc    6%
  fx_advance            2.22 ms   1.33 Mcyc    3%
  fx_timeline_step     63.5 us   38.13 kcyc    0%
  canvas_clear         88.9 us   53.35 kcyc    0%
  canvas_buffer_wait   11.00 ms   6.60 Mcyc   17%
```

Wall min/avg/max = 60.38/62.45/64.43 ms. The composed shader draw owns 72%
of wall time; preparation and runtime advance add 6.09 ms. Display-sync idle
absorbs the remaining 11.00 ms on average.

### Per-pixel figures

The draw covers one 10,368-sample quadrant. The worst window spends about
2,616 cycles per sample in `fx_shader_draw`; no per-pixel profiling scope is
enabled, so that figure includes the whole composed pipeline.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.78/2.00/17.17 us  cpu 3.68%
isr_pack         144/frame  min/avg/max 6.32/7.10/9.28 us   cpu 1.63%
isr_dma_submit   144/frame  min/avg/max 0.67/0.94/2.99 us   cpu 0.21%
```

- Packing dominates submit CPU cost; wire transfer remains asynchronous.
- The measured ISR share is included in the free-running cycle scopes.
- The 52.39 ms render peak stays below the 62.5 ms cadence budget.

## Summary ranking

1. `fx_shader_draw` — 72% of the worst window, 45.20 ms/frame.
2. `fx_prepare_frame` — 6%, 3.87 ms/frame.
3. `fx_advance` — 3%, 2.22 ms/frame.

The global-O3 twin peaks at 48.35 ms, 7.7% below the shipping peak.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- No per-pixel scope was enabled, avoiding instrumentation distortion.
- The shipping image uses the repository's selective `HS_O3` regions on an
  otherwise `-Os` build.
- The source tree had no tracked WIP; review/profile documents are ignored
  working-tree artifacts.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AshCloud`,
`HS_PROFILE_WINDOW=32`; `just profile AshCloud` builds, flashes and captures.
