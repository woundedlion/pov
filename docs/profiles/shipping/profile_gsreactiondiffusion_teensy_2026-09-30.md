# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-09-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).

Raw capture: [preserved capture](../evidence/gs_optimization_2026-09-30/finalship/capture.txt), originally `build/prof/gsreactiondiffusion_ship.log`, captured 2026-09-30 01:19 on COM3. Replaces `profile_gsreactiondiffusion_teensy_2026-09-28.md`.

[Optimization campaign](../gsreactiondiffusion_optimization_2026-09-30.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base plus selective `HS_O3_FN` physics, shader, cull, lattice-distance/refinement and orientation kernels |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, source tip `0fa3bede8ca9d4f6ecb6afc94deb10429d29e087` (clean captured source tree) |
| Method | `HS_PROFILE`, 32-frame windows, 130 s capture, epoch stretched to 1200 revolutions (150 s). Exact runtime frames 2–2048; frame 1 excluded. Scope/ISR summaries use complete windows 33–2048. |
| Reproduce | `HS_TEENSY_PORT=COM3 HS_PROFILE_TREE=/c/work/Holosphere bash tools/profile_one.sh GSReactionDiffusion profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Single-effect image size: `FLASH: code:60936, data:251868, headers:8732   free for files:1710080` / `RAM1: variables:315168, code:23800, padding:8968   free for local variables:176352` / `RAM2: variables:520064  free for malloc/new:4224`.

Full-roster shipping Phantasm size, separately built with profiling disabled: `FLASH: code:529648, data:746232, headers:8216   free for files:747520` / `RAM1: variables:314784, code:166680, padding:29928   free for local variables:12896` / `RAM2: variables:520064  free for malloc/new:4224`. This is the shipping RAM1/ITCM budget; the single-effect image above is not the roster memory budget.

Exactness cross-check: window frames 449–480, root 1,198,957,710 cycles ÷ 600 MHz versus wall sum 1,998,267 us, within **2.1 ppm**. [Parser validation](../evidence/gs_optimization_2026-09-30/finalship/validate.txt) reports VALID.

## Frame cadence

**Runtime aggregate** (frame 1 excluded): mean render **32.241 ms/frame**, peak **36.528 ms** at frame 464, spilled **0/2047 (0.00%)**. Runtime mean wall is 62.432 ms/frame. Scope totals below cover frames 33–2048 only.

Startup setup render: **23.542 ms** at frame 1, before publication and excluded from runtime statistics.

A display window is 62.5 ms. The effect renders one 144×72 quadrant (10,368 pixels), with four samples per shaded pixel. Peak render leaves 25.972 ms against this deadline. Every captured live frame holds the 16 fps render budget. `canvas_buffer_wait` is idle time until the next display flip.

| Scope | Mean ms/frame | Worst window ms/frame | Frames |
|---|--:|--:|---|
| `grd_render` | 32.154 | 36.050 | 449–480 |
| `grd_simulate` | 6.712 | 6.716 | 1889–1920 |
| `grd_rasterize` | 25.330 | 28.797 | 449–480 |
| `grd_shader_draw` | 24.448 | 28.103 | 449–480 |
| `grd_cull_flags` | 0.507 | 0.952 | 1505–1536 |
| `grd_orient` | 0.375 | 0.376 | 161–192 |
| `canvas_buffer_wait` | 30.041 | 44.557 | 993–1024 |

## Phase-by-phase readout

Phase schedule: reaction growth changes the active area; settling triggers dissolve and reseeding. This capture has no lifecycle markers, so windows below describe measured workload rather than assigning unobserved lifecycle boundaries.

### Early captured growth (frames 33–64)

```
frame                      62.65 ms  37.59 Mcyc 100.0%
  pov_preserve_half        139.7 us   83.8 kcyc   0.2% x1 139.7us/call
  grd_render               27.90 ms  16.74 Mcyc  44.5%
    grd_rasterize          21.19 ms  12.71 Mcyc  33.8%
      grd_shader_draw      20.03 ms  12.02 Mcyc  32.0% x1 20025.7us/call
      grd_cull_flags       784.0 us  470.4 kcyc   1.3% x1 784.0us/call
      grd_orient           375.5 us  225.3 kcyc   0.6% x1 375.5us/call
    grd_simulate            6.71 ms   4.03 Mcyc  10.7% x1 6710.0us/call
  rd_timeline_step          38.0 us   22.8 kcyc   0.1% x1 38.0us/call
  canvas_clear              84.5 us   50.7 kcyc   0.1% x1 84.5us/call
  canvas_buffer_wait       34.49 ms  20.69 Mcyc  55.1% x1 34490.0us/call
```

Wall min/avg/max = 57.951/62.649/67.149 ms. `grd_shader_draw` costs 20.026 ms/frame and simulation 6.710 ms/frame. Increased active coverage raises shader work; display synchronization absorbs the remaining time. These scope figures are window averages, not individual-frame peaks.

### Highest mean render window (frames 449–480)

```
frame                      62.45 ms  37.47 Mcyc 100.0%
  pov_preserve_half        138.5 us   83.1 kcyc   0.2% x1 138.5us/call
  grd_render               36.05 ms  21.63 Mcyc  57.7%
    grd_rasterize          28.80 ms  17.28 Mcyc  46.1%
      grd_shader_draw      28.10 ms  16.86 Mcyc  45.0% x1 28102.9us/call
      grd_cull_flags       323.1 us  193.9 kcyc   0.5% x1 323.1us/call
      grd_orient           370.8 us  222.5 kcyc   0.6% x1 370.8us/call
    grd_simulate            6.72 ms   4.03 Mcyc  10.8% x1 6715.8us/call
  rd_timeline_step          29.6 us   17.8 kcyc   0.0% x1 29.6us/call
  canvas_clear              84.3 us   50.6 kcyc   0.1% x1 84.3us/call
  canvas_buffer_wait       26.14 ms  15.69 Mcyc  41.9% x1 26141.8us/call
```

Wall min/avg/max = 62.057/62.446/62.814 ms. `grd_shader_draw` costs 28.103 ms/frame and simulation 6.716 ms/frame. Increased active coverage raises shader work; display synchronization absorbs the remaining time. These scope figures are window averages, not individual-frame peaks.

### Per-pixel figures

No `filter_blend` counter is present. The raster covers 10,368 quadrant positions, but hot-pixel counts are unavailable in this standard capture; dividing by the full quadrant would not measure cost per shaded pixel.

## Column-ISR / DMA marshaling cost

Complete windows 33–2048; columns are calls/frame, per-call min/avg/max, and CPU share. Pack and submit are nested inside wake and must not be added to its inclusive CPU share.

```
isr_wake          1152.1/f 0.61/1.69/17.02 us  3.12% CPU
  isr_pack         144.0/f 6.24/7.06/9.94 us  1.63% CPU
  isr_dma_submit   144.0/f 0.60/0.93/9.68 us  0.21% CPU
```

- Pack averages 7.062 us/call; DMA submission averages 0.930 us/call. CPU shares use capture-window elapsed time; averages use summed logged ISR microseconds/counts.
- SPI DMA proceeds asynchronously. At 24 MHz, the 72-LED 300-byte frame takes approximately 115 us including configured LPSPI framing; an image plus trailing black frame takes 230 us. These are transport calculations, not measured CPU submission times.
- Inclusive wake share 3.12% leaves approximately 60.552 ms of foreground CPU per 62.5 ms window before other interrupts. Measured render already includes interrupt time; no further speedup is required to meet the observed deadline.

## Summary ranking

1. `grd_shader_draw` — 45.0% of the highest-render window, 28.103 ms/frame.
2. `grd_simulate` — 10.8% of the highest-render window, 6.716 ms/frame.
3. `grd_orient` — 0.6% of the highest-render window, 0.371 ms/frame.
4. `grd_cull_flags` — 0.5% of the highest-render window, 0.323 ms/frame.

Native runs validate numerical behavior; no directly comparable WASM/native timing ledger entry is available for this exact revision.

## Caveats

- All foreground scopes include ISR time because CYCCNT free-runs. Wake includes nested pack/submit costs.
- `filter_blend` normally parents under whichever scope first enters it and can disappear with an inactive parent; this capture has no such counter.
- No per-pixel profiling scopes are enabled. The separate diagnostic build adds overhead and is not used for these timings.
- Shipping selective-O3 covers GS `step_physics`, `shade_pixel`, `fill_hot_flags`, shared lattice distance/refinement and orientation methods, plus the renderer/driver hot paths. Global-O3 changes the single-effect comparison image globally; it is not the full-roster shipping configuration.
- The epoch is extended to avoid reinitialization during capture; there is no dwell compression or simulation-speed override.
- Captured from a clean source tree; the archived source diff is empty. The raw provenance footer and [source patch](../evidence/gs_optimization_2026-09-30/finalship/source.json) identify the measured working tree; the evidence retains build sizes and hashes.
- Frame 1 is excluded from exact runtime rows, while the entire first window is excluded from scope and ISR aggregates. ISR totals are logged in integer microseconds, so their aggregate per-call averages have small quantization error.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`. `just profile GSReactionDiffusion` is the basic shortcut; use the Setup command to reproduce this duration and epoch under the shared-device lock.
