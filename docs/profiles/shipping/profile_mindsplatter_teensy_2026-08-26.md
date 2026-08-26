# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `build/prof/mindsplatter_ship.log`, sourced from the isolated
regression-reclaim tree. Replaces the earlier same-day sweep capture after
the cadence-reclamation changes in `0df961b81`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses direct-AA sink, Plot cull/raster, and splat `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, tip `0df961b818ae08c5f58639edb93d812164a9356a` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 110 s capture |
| Reproduce | `bash tools/profile_one.sh MindSplatter profile 110 16` |

Image size:

```text
FLASH: code:65,616, data:545,256, headers:8,648
       free for files:1,412,096
RAM1:  variables:315,264, code:42,248, padding:23,288
       free for local variables:143,488
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 993–1008 root counter cycles ÷ 600 MHz
match the measured wall sum within **3.7 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, the closed preset cycle, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`msp_draw_particles` averages 20.70 ms/frame; its worst window is
46.30 ms/frame (frames 993–1008). Peak frame render is **52.77 ms**
(frames 993–1008); spilled **0/1728 frames (0.0%)**.

A display window is 62.5 ms. The segmented master renders one quadrant,
approximately 10,368 pixels. Every captured preset holds 16 fps; the peak
retains 9.73 ms of the display window. The `canvas_buffer_wait` scope is
round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: the parser emits nine ownership buckets (the initial marker
state plus the eight presets), confirms all 8/8 presets, and closes back to
preset 0. The sections below bound the clean-hold regimes; preset ownership,
including the following transition, is tabulated separately.

### Highest-cost clean hold (window frames 993–1008)

```text
frame                        62.18ms 37.31Mcyc 100%
  msp_draw_particles         46.30ms 27.78Mcyc  74%
    msp_particle_scan        46.29ms 27.77Mcyc  74%
      plot_ps_raster         32.80ms 19.68Mcyc  53% x557.8 58.8us/c
      plot_ps_deferred       559.6us 335.80kcyc   1% x557.8 1.0us/c
      plot_ps_gate            7.11ms  4.26Mcyc  11%
        plot_ps_cartesian_gate 1.06ms 633.68kcyc 2% x1537 0.7us/c
      plot_ps_tween           4.98ms  2.99Mcyc   8% x1537 3.2us/c
  msp_particle_step           3.10ms  1.86Mcyc   5% x1.0 3097.1us/c
  msp_timeline_step           68.1us 40.87kcyc   0% x1.0 68.1us/c
  canvas_clear                89.2us 53.52kcyc   0% x1.0 89.2us/c
  canvas_buffer_wait         12.63ms  7.58Mcyc  20% x1.0 12625.9us/c
```

Wall min/avg/max = 59.46/62.18/64.39 ms. `msp_draw_particles` accounts for
74.5% of this measured frame at 46.30 ms/frame. Complete render averages
49.55 ms; `canvas_buffer_wait` contributes 12.63 ms of flip-alignment idle.
The preset's exact peak render remains below one display window.

### Startup / lowest-cost regime (window frames 1–16)

```text
frame                        54.75ms 32.85Mcyc 100%
  msp_draw_particles         419.9us 251.94kcyc   1%
    msp_particle_scan        412.6us 247.59kcyc   1%
      plot_ps_raster         206.7us 124.02kcyc   0% x15.0 13.8us/c
      plot_ps_deferred         5.6us  3.41kcyc    0% x15.0 0.4us/c
      plot_ps_gate            87.8us 52.71kcyc    0%
        plot_ps_cartesian_gate 21.9us 13.16kcyc 0% x60.0 0.4us/c
      plot_ps_tween           78.2us 46.93kcyc    0% x60.0 1.3us/c
  msp_particle_step          132.3us 79.40kcyc    0% x1.0 132.3us/c
  msp_timeline_step           41.6us 24.96kcyc    0% x1.0 41.6us/c
  canvas_clear                86.9us 52.20kcyc    0% x1.0 86.9us/c
  canvas_buffer_wait         54.07ms 32.44Mcyc   99% x1.0 54067.3us/c
```

Wall min/avg/max = 0.19/54.75/62.64 ms. Startup has not yet populated the
particle field: draw time is 0.42 ms/frame and complete render averages
0.68 ms. The 54.07 ms wait is display-sync idle, not rendering work.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to preset 0.
Rows rank the modal-call-count clean holds by `msp_draw_particles` cost.
Cadence peaks use per-frame ownership, including the transition following each
preset. All 108 owned windows meet the modal clean-hold criterion.

| rank | preset | windows | blended px/f | `msp_draw_particles` ms | render ms | fps |
|---:|---:|--:|--:|--:|--:|--:|
| 1 | `7` | 10/10 | — | 46.30 | 52.77 | 16.0 |
| 2 | `6` | 10/10 | — | 40.08 | 50.36 | 16.0 |
| 3 | `8` | 10/10 | — | 28.89 | 42.30 | 16.0 |
| 4 | `2` | 20/20 | — | 28.30 | 34.16 | 16.0 |
| 5 | `3` | 18/18 | — | 27.26 | 41.36 | 16.0 |
| 6 | `1` | 10/10 | — | 27.21 | 32.34 | 16.0 |
| 7 | `4` | 10/10 | — | 26.94 | 31.68 | 16.0 |
| 8 | `0` | 10/10 | — | 19.59 | 22.97 | 16.0 |
| 9 | `5` | 10/10 | — | 17.10 | 21.61 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the highest-cost window; the path
writes directly. `msp_draw_particles` uses 2679.1 cycles per configured
quadrant sample over 10,368 samples/frame; `plot_ps_raster` accounts for
1898.1 of those cycles/sample.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1156.3/f  0.54/1.65/22.59us  cpu 3.06%
isr_pack          143.9/f  6.23/6.88/9.66us   cpu 1.59%
isr_dma_submit    143.9/f  0.59/0.94/1.97us   cpu 0.22%
```

- Pack plus submit costs 1.12 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling
  only.
- Aggregate ISR CPU share is 4.86% of one 62.5 ms window. The peak render
  requires 0.84× the one-window budget, so the cadence target needs no
  speedup.

## Summary ranking

1. `msp_draw_particles` — 33.2% of aggregate root time, 20.70 ms/frame.
2. `msp_particle_scan` — 33.2% of aggregate root time, 20.69 ms/frame.
3. `plot_ps_raster` — 21.2% of aggregate root time, 13.20 ms/frame.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches direct-AA sink, Plot cull/raster, and splat
  `HS_O3` regions; the rest of the image retains the `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each
  preset, not its per-frame cost. No capture-time compression knob was used.
- Presets 2/8, 3/8, and 7/8 select Cube geometry. Emitters remain uncapped and
  every active emitter continues to spawn each frame.
- Provenance attests a clean source tree at
  `0df961b818ae08c5f58639edb93d812164a9356a`; no uncommitted source state was
  profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=MindSplatter` and
`HS_PROFILE_WINDOW=16`; `just profile MindSplatter` performs the locked build,
flash, capture, and artifact attestation.
