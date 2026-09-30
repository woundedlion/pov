# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-09-30, **-O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).

[Shipping selective-O3 sibling](../shipping/profile_gsreactiondiffusion_teensy_2026-09-30.md).

Raw capture: [preserved capture](../evidence/gs_optimization_2026-09-30/finalo3/capture.txt), originally `build/prof/gsreactiondiffusion_o3.log`, captured 2026-09-30 01:23 on COM3. Replaces `profile_gsreactiondiffusion_teensy_2026-08-26.md`.

[Optimization campaign](../gsreactiondiffusion_optimization_2026-09-30.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3` env: `-O3 -ffast-math` globally, single-effect comparison image |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, source tip `3968cd8f65f7d21ad6f56597116a0cf49c13a17c` (clean captured source tree) |
| Method | `HS_PROFILE`, 32-frame windows, 130 s capture, epoch stretched to 1200 revolutions (150 s). Exact runtime frames 2–2048; frame 1 excluded. Scope/ISR summaries use complete windows 33–2048. |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=/c/work/Holosphere bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Single-effect image size: `FLASH: code:74792, data:251864, headers:8192   free for files:1696768` / `RAM1: variables:315200, code:32696, padding:72   free for local variables:176320` / `RAM2: variables:520064  free for malloc/new:4224`.

Full-roster shipping Phantasm size, separately built with profiling disabled: `FLASH: code:529648, data:746232, headers:8216   free for files:747520` / `RAM1: variables:314784, code:166680, padding:29928   free for local variables:12896` / `RAM2: variables:520064  free for malloc/new:4224`. This is the shipping RAM1/ITCM budget; the single-effect image above is not the roster memory budget.

Exactness cross-check: window frames 449–480, root 1,198,913,385 cycles ÷ 600 MHz versus wall sum 1,998,190 us, within **0.5 ppm**. [Parser validation](../evidence/gs_optimization_2026-09-30/finalo3/validate.txt) reports VALID.

## Frame cadence

**Runtime aggregate** (frame 1 excluded): mean render **31.873 ms/frame**, peak **36.477 ms** at frame 464, spilled **0/2047 (0.00%)**. Runtime mean wall is 62.427 ms/frame. Scope totals below cover frames 33–2048 only.

Startup setup render: **23.615 ms** at frame 1, before publication and excluded from runtime statistics.

A display window is 62.5 ms. The effect renders one 144×72 quadrant (10,368 pixels), with four samples per shaded pixel. Peak render leaves 26.023 ms against this deadline. Every captured live frame holds the 16 fps render budget. `canvas_buffer_wait` is idle time until the next display flip.

| Scope | Mean ms/frame | Worst window ms/frame | Frames |
|---|--:|--:|---|
| `grd_render` | 31.782 | 35.952 | 449–480 |
| `grd_simulate` | 6.630 | 6.634 | 1889–1920 |
| `grd_rasterize` | 25.044 | 28.794 | 449–480 |
| `grd_shader_draw` | 24.161 | 28.096 | 449–480 |
| `grd_cull_flags` | 0.509 | 0.949 | 1505–1536 |
| `grd_orient` | 0.373 | 0.374 | 385–416 |
| `canvas_buffer_wait` | 30.409 | 42.595 | 993–1024 |

## Phase-by-phase readout

Phase schedule: reaction growth changes the active area; settling triggers dissolve and reseeding. This capture has no lifecycle markers, so windows below describe measured workload rather than assigning unobserved lifecycle boundaries.

### Early captured growth (frames 33–64)

```
frame                      62.64 ms  37.59 Mcyc 100.0%
  pov_preserve_half        141.6 us   85.0 kcyc   0.2% x1 141.6us/call
  grd_render               27.83 ms  16.70 Mcyc  44.4%
    grd_rasterize          21.20 ms  12.72 Mcyc  33.8%
      grd_shader_draw      20.04 ms  12.02 Mcyc  32.0% x1 20037.9us/call
      grd_cull_flags       787.6 us  472.5 kcyc   1.3% x1 787.6us/call
      grd_orient           373.5 us  224.1 kcyc   0.6% x1 373.5us/call
    grd_simulate            6.63 ms   3.98 Mcyc  10.6% x1 6626.9us/call
  rd_timeline_step          37.7 us   22.6 kcyc   0.1% x1 37.7us/call
  canvas_clear              84.2 us   50.5 kcyc   0.1% x1 84.2us/call
  canvas_buffer_wait       34.55 ms  20.73 Mcyc  55.2% x1 34551.2us/call
```

Wall min/avg/max = 58.009/62.642/67.100 ms. `grd_shader_draw` costs 20.038 ms/frame and simulation 6.627 ms/frame. Increased active coverage raises shader work; display synchronization absorbs the remaining time. These scope figures are window averages, not individual-frame peaks.

### Highest mean render window (frames 449–480)

```
frame                      62.44 ms  37.47 Mcyc 100.0%
  pov_preserve_half        139.3 us   83.6 kcyc   0.2% x1 139.3us/call
  grd_render               35.95 ms  21.57 Mcyc  57.6%
    grd_rasterize          28.79 ms  17.28 Mcyc  46.1%
      grd_shader_draw      28.10 ms  16.86 Mcyc  45.0% x1 28096.1us/call
      grd_cull_flags       326.3 us  195.8 kcyc   0.5% x1 326.3us/call
      grd_orient           371.8 us  223.1 kcyc   0.6% x1 371.8us/call
    grd_simulate            6.63 ms   3.98 Mcyc  10.6% x1 6634.1us/call
  rd_timeline_step          28.7 us   17.2 kcyc   0.0% x1 28.7us/call
  canvas_clear              84.4 us   50.6 kcyc   0.1% x1 84.4us/call
  canvas_buffer_wait       26.24 ms  15.74 Mcyc  42.0% x1 26238.0us/call
```

Wall min/avg/max = 62.047/62.443/62.914 ms. `grd_shader_draw` costs 28.096 ms/frame and simulation 6.634 ms/frame. Increased active coverage raises shader work; display synchronization absorbs the remaining time. These scope figures are window averages, not individual-frame peaks.

### Per-pixel figures

No `filter_blend` counter is present. The raster covers 10,368 quadrant positions, but hot-pixel counts are unavailable in this standard capture; dividing by the full quadrant would not measure cost per shaded pixel.

## Column-ISR / DMA marshaling cost

Complete windows 33–2048; columns are calls/frame, per-call min/avg/max, and CPU share. Pack and submit are nested inside wake and must not be added to its inclusive CPU share.

```
isr_wake          1152.1/f 0.33/1.52/15.82 us  2.80% CPU
  isr_pack         144.0/f 6.00/6.93/9.78 us  1.60% CPU
  isr_dma_submit   144.0/f 0.60/0.93/6.42 us  0.21% CPU
```

- Pack averages 6.929 us/call; DMA submission averages 0.929 us/call. CPU shares use capture-window elapsed time; averages use summed logged ISR microseconds/counts.
- SPI DMA proceeds asynchronously. At 24 MHz, the 72-LED 300-byte frame takes approximately 115 us including configured LPSPI framing; an image plus trailing black frame takes 230 us. These are transport calculations, not measured CPU submission times.
- Inclusive wake share 2.80% leaves approximately 60.748 ms of foreground CPU per 62.5 ms window before other interrupts. Measured render already includes interrupt time; no further speedup is required to meet the observed deadline.

## Summary ranking

1. `grd_shader_draw` — 45.0% of the highest-render window, 28.096 ms/frame.
2. `grd_simulate` — 10.6% of the highest-render window, 6.634 ms/frame.
3. `grd_orient` — 0.6% of the highest-render window, 0.372 ms/frame.
4. `grd_cull_flags` — 0.5% of the highest-render window, 0.326 ms/frame.

Native runs validate numerical behavior; no directly comparable WASM/native timing ledger entry is available for this exact revision.

## Caveats

- All foreground scopes include ISR time because CYCCNT free-runs. Wake includes nested pack/submit costs.
- `filter_blend` normally parents under whichever scope first enters it and can disappear with an inactive parent; this capture has no such counter.
- No per-pixel profiling scopes are enabled. The separate diagnostic build adds overhead and is not used for these timings.
- Shipping selective-O3 covers GS `step_physics`, `shade_pixel`, `fill_hot_flags`, shared lattice distance/refinement and orientation methods, plus the renderer/driver hot paths. Global-O3 changes the single-effect comparison image globally; it is not the full-roster shipping configuration.
- The epoch is extended to avoid reinitialization during capture; there is no dwell compression or simulation-speed override.
- Captured from a clean source tree; the archived source diff is empty. The raw provenance footer and [source patch](../evidence/gs_optimization_2026-09-30/finalo3/source.json) identify the measured working tree; the evidence retains build sizes and hashes.
- Frame 1 is excluded from exact runtime rows, while the entire first window is excluded from scope and ISR aggregates. ISR totals are logged in integer microseconds, so their aggregate per-call averages have small quantization error.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`. `just profile GSReactionDiffusion` is the basic shortcut; use the Setup command to reproduce this duration and epoch under the shared-device lock.

### Global -O3 vs selective -O3

Shipping peak 36.528 ms versus global-O3 36.477 ms: 1.001× peak speed ratio. Single-effect global-O3 minus shipping image: FLASH code **+13,856 B**, ITCM **+8,896 B**. These deltas compare the two single-effect profile images, not Phantasm against a profile image.
