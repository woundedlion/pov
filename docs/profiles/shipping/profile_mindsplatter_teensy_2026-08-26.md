# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `build/prof/mindsplatter_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_mindsplatter_teensy_2026-08-07.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses direct-AA sink, Plot cull/raster, and splat `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile 110 16` |

Image size:

```text
FLASH: code:65,656, data:545,256, headers:8,608
       free for files:1,412,096
RAM1:  variables:315,264, code:42,280, padding:23,256
       free for local variables:143,488
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 1633–1648 root counter
cycles ÷ 600 MHz match the measured wall sum within **1.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `msp_draw_particles` averages
24.49 ms/frame; its worst window is 51.73 ms/frame
(frames 1633–1648). Peak frame render is
**67.05 ms** (frames 1633–1648);
spilled **39/1696 frames**
(2.3%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak exceeds one display window by 4.55 ms. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 1633–1648)

```text
frame                        108.90 ms  65.34 Mcyc  100%
  msp_draw_particles          51.73 ms  31.04 Mcyc   48%
    msp_particle_scan         51.73 ms  31.04 Mcyc   47%
      plot_ps_raster          36.48 ms  21.89 Mcyc   33%  x520.1 70.1us/c
      plot_ps_deferred         3.28 ms   1.97 Mcyc    3%  x520.1 6.3us/c
      plot_ps_gate             6.50 ms   3.90 Mcyc    6%
        plot_ps_cartesian_g   985.4 us 591.26 kcyc    1%  x1475.9 0.7us/c
      plot_ps_tween            4.67 ms   2.80 Mcyc    4%  x1475.9 3.2us/c
  msp_particle_step           11.25 ms   6.75 Mcyc   10%  x1.0 11246.1us/c
  msp_timeline_step            62.0 us  37.21 kcyc    0%  x1.0 62.0us/c
  canvas_clear                 85.8 us  51.51 kcyc    0%  x1.0 85.8us/c
  canvas_buffer_wait          45.77 ms  27.46 Mcyc   42%  x1.0 45772.1us/c
```

Wall min/avg/max = 55.25/108.90/126.14 ms. `msp_draw_particles`
accounts for 47.5% of this measured frame at 51.73 ms/frame.
The complete render is 63.13 ms; `canvas_buffer_wait` contributes
45.77 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–16)

```text
frame                         54.75 ms  32.85 Mcyc  100%
  msp_draw_particles          426.8 us 256.06 kcyc    1%
    msp_particle_scan         419.5 us 251.71 kcyc    1%
      plot_ps_raster          212.2 us 127.34 kcyc    0%  x15.0 14.1us/c
      plot_ps_deferred          5.7 us   3.43 kcyc    0%  x15.0 0.4us/c
      plot_ps_gate             90.0 us  54.01 kcyc    0%
        plot_ps_cartesian_g    22.4 us  13.49 kcyc    0%  x60.0 0.4us/c
      plot_ps_tween            78.6 us  47.15 kcyc    0%  x60.0 1.3us/c
  msp_particle_step           131.0 us  78.61 kcyc    0%  x1.0 131.0us/c
  msp_timeline_step            40.4 us  24.28 kcyc    0%  x1.0 40.4us/c
  canvas_clear                 86.8 us  52.05 kcyc    0%  x1.0 86.8us/c
  canvas_buffer_wait          54.06 ms  32.44 Mcyc   99%  x1.0 54060.7us/c
```

Wall min/avg/max = 0.18/54.75/62.63 ms. `msp_draw_particles`
accounts for 0.8% of this measured frame at 0.43 ms/frame.
The complete render is 0.69 ms; `canvas_buffer_wait` contributes
54.06 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `msp_draw_particles` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `msp_draw_particles` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `3` | — | 16/16 | — | 51.73 | 67.05 | 8.0 |
| 2 | `7` | — | 10/10 | — | 50.48 | 64.61 | 8.0 |
| 3 | `2` | — | 20/20 | — | 43.74 | 60.54 | 16.0 |
| 4 | `6` | — | 10/10 | — | 40.73 | 51.68 | 16.0 |
| 5 | `4` | — | 10/10 | — | 33.13 | 39.78 | 16.0 |
| 6 | `1` | — | 10/10 | — | 28.50 | 33.82 | 16.0 |
| 7 | `8` | — | 10/10 | — | 27.89 | 39.26 | 16.0 |
| 8 | `0` | — | 10/10 | — | 20.05 | 23.45 | 16.0 |
| 9 | `5` | — | 10/10 | — | 17.31 | 23.11 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `msp_draw_particles` uses 2993.9 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1182.3/f  0.48/1.67/20.48 us  cpu 3.17%
isr_pack          147.1/f  6.23/6.91/9.66 us  cpu 1.63%
isr_dma_submit    147.1/f  0.60/0.94/2.31 us  cpu 0.22%
```

- Pack plus submit costs 1.15 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.01% of one 62.5 ms window. The peak
  render requires 1.07× the one-window budget, so the cadence target
  needs a 1.07× speedup.

## Summary ranking

1. `msp_draw_particles` — 38.4% of aggregated root time, 24.49 ms/frame: inclusive measured scope in the live driver.
2. `msp_particle_scan` — 38.4% of aggregated root time, 24.49 ms/frame: inclusive measured scope in the live driver.
3. `plot_ps_raster` — 24.3% of aggregated root time, 15.52 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches direct-AA sink, Plot cull/raster, and splat `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=MindSplatter` and
`HS_PROFILE_WINDOW=16`; `just profile MindSplatter` performs the
locked build, flash, capture, and artifact attestation.
