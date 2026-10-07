# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-09-30, **-O3**)

Point-in-time snapshot (regenerate with the explicit `profile_o3` Reproduce command).

The same-source shipping sibling has been superseded by the
[current shipping capture](../shipping/profile_gsreactiondiffusion_teensy_2026-10-06.md).
This O3 report retains its September 30 source and image comparisons.

Raw capture: preserved capture (supporting artifact removed), captured 2026-09-30 19:15 on COM4. Replaces the earlier 2026-09-30 report from before pigment blending, hue rotation and shimmer were enabled. Optimization campaign and current-state update (supporting artifact removed).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: `-O3 -ffast-math` globally, single-effect comparison image |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, source `7baf3cc430753e1a3693134150b282f212aa4e9e` (clean captured candidate) |
| Method | `HS_PROFILE`, 32-frame windows, 130 s capture, epoch stretched to 1200 revolutions (150 s). Exact runtime frames 2–539; setup frame 1 excluded. Scope/ISR summaries use complete windows 33–512. |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE=/c/work/Holosphere bash tools/profile_one.sh GSReactionDiffusion profile_o3 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Single-effect image size: `FLASH: code:82376, data:252244, headers:8420   free for files:1688576` / `RAM1: variables:315136, code:38872, padding:26664   free for local variables:143616` / `RAM2: variables:520064  free for malloc/new:4224`.

Full-roster shipping Phantasm, separately built with profiling disabled: `FLASH: code:535552, data:746852, headers:8860   free for files:740352` / `RAM1: variables:314784, code:171928, padding:24680   free for local variables:12896` / `RAM2: variables:520064  free for malloc/new:4224`. The full-roster build passes the region-budget and layout gates; its RAM1 code is the shipping ITCM budget, separate from the single-effect comparison image.

Exactness cross-check: window frames 385–416, root 5,833,189,292 cycles ÷ 600 MHz versus measured wall sum 9,721,985 us, within **0.3 ppm**. Parser validation (supporting artifact removed) reports VALID.

## Frame cadence

README cells: peak 🔴 279.686, spilled 🔴 538/538 (100.00%).

**Runtime aggregate** (setup frame 1 excluded): mean render **212.279 ms/frame**, peak **279.686 ms** at frame 395, spilled **538/538 (100.00%)**. Runtime mean wall is 238.277 ms/frame, or 4.20 rendered frames/s across these live rows.

Startup setup render: **127.055 ms** at frame 1, before publication and excluded from runtime statistics.

A display window is 62.5 ms. The effect renders one 144×72 quadrant (10,368 positions), with four samples per shaded pixel. Every captured live frame exceeds the 16 fps render deadline. Peak render exceeds that deadline by 217.186 ms and requires 4.47× render speedup to fit one window. `canvas_buffer_wait` is synchronization idle until a display flip, not rendering work.

| Scope | Mean ms/frame | Worst window ms/frame |
|---|--:|--:|
| `grd_render` | 218.408 | 266.581 |
| `grd_simulate` | 93.779 | 93.871 |
| `grd_rasterize` | 121.494 | 169.294 |
| `grd_shader_draw` | 120.561 | 168.559 |
| `grd_cull_flags` | 0.533 | 0.956 |
| `grd_orient` | 0.400 | 0.401 |
| `canvas_buffer_wait` | 24.608 | 46.022 |

## Phase-by-phase readout

Phase schedule: the evolving reaction field changes active shading coverage. The capture has no lifecycle markers, so these windows describe measured workloads without assigning unobserved dissolve/reseed boundaries. The slower image advances fewer simulation frames within 130 seconds; this is a measured interval, not a lifetime worst-case bound.

### Early captured growth (frames 33–64)

```
frame                     188.55 ms 113.13 Mcyc 100.0%
  pov_preserve_half        133.8 us   80.3 kcyc   0.1% x1 133.8us/call
  grd_render              157.73 ms  94.64 Mcyc  83.7%
    grd_rasterize          61.04 ms  36.63 Mcyc  32.4%
      grd_shader_draw      59.85 ms  35.91 Mcyc  31.7% x1 59848.7us/call
      grd_cull_flags       798.8 us  479.3 kcyc   0.4% x1 798.8us/call
      grd_orient           396.4 us  237.8 kcyc   0.2% x1 396.4us/call
    grd_simulate           93.70 ms  56.22 Mcyc  49.7% x1 93700.5us/call
  rd_timeline_step          37.3 us   22.4 kcyc   0.0% x1 37.3us/call
  canvas_clear              84.3 us   50.6 kcyc   0.0% x1 84.3us/call
  canvas_buffer_wait       30.56 ms  18.34 Mcyc  16.2% x1 30563.5us/call
```

Wall min/avg/max = 150.422/188.548/226.606 ms. Shader work costs 59.849 ms/frame and simulation 93.700 ms/frame. Both contribute materially to the deadline overrun; higher active coverage raises the shader cost. These scope costs are window means, distinct from the exact individual-frame peak.

### Highest mean render window (frames 385–416)

```
frame                     303.81 ms 182.29 Mcyc 100.0%
  pov_preserve_half        131.2 us   78.7 kcyc   0.0% x1 131.2us/call
  grd_render              266.58 ms 159.95 Mcyc  87.7%
    grd_rasterize         169.29 ms 101.58 Mcyc  55.7%
      grd_shader_draw     168.56 ms 101.14 Mcyc  55.5% x1 168558.7us/call
      grd_cull_flags       336.7 us  202.0 kcyc   0.1% x1 336.7us/call
      grd_orient           398.6 us  239.2 kcyc   0.1% x1 398.6us/call
    grd_simulate           93.77 ms  56.26 Mcyc  30.9% x1 93767.9us/call
  rd_timeline_step          31.3 us   18.8 kcyc   0.0% x1 31.3us/call
  canvas_clear              85.0 us   51.0 kcyc   0.0% x1 85.0us/call
  canvas_buffer_wait       36.98 ms  22.19 Mcyc  12.2% x1 36983.2us/call
```

Wall min/avg/max = 246.890/303.812/332.039 ms. Shader work costs 168.559 ms/frame and simulation 93.768 ms/frame. Both contribute materially to the deadline overrun; higher active coverage raises the shader cost. These scope costs are window means, distinct from the exact individual-frame peak.

### Per-pixel figures

No `filter_blend` counter is present. The raster visits 10,368 quadrant positions, but active-pixel counts are unavailable. Dividing by the full quadrant would not measure cost per shaded pixel; no per-pixel cost is inferred.

## Column-ISR / DMA marshaling cost

Complete post-startup windows only; columns are calls/frame, per-call min/avg/max, and CPU share. Pack and submit are nested inside wake, so their shares are not added to its inclusive share.

```
isr_wake          4485.1/f 0.41/1.59/20.34 us  2.94% CPU
  isr_pack         560.6/f 5.98/6.94/11.06 us  1.60% CPU
  isr_dma_submit   560.6/f 0.60/0.93/10.93 us  0.22% CPU
```

- Pack averages 6.943 us/call; DMA submission averages 0.935 us/call. CPU shares use logged window elapsed time.
- SPI DMA transfers asynchronously. At 24 MHz, a 72-LED, 300-byte frame takes about 115 us including LPSPI framing; image plus trailing black frame takes about 230 us. These are wire-time calculations, not measured CPU submission times.
- Inclusive wake share 2.94% leaves about 60.666 ms of foreground CPU per 62.5 ms interval before other interrupts. Foreground scope times already include interrupts; the observed peak requires 4.47× render speedup to fit the wall deadline.

## Summary ranking

1. `grd_shader_draw` — 55.5% of the highest-render window, 168.559 ms/frame.
2. `grd_simulate` — 30.9% of the highest-render window, 93.768 ms/frame.
3. `grd_orient` — 0.1% of the highest-render window, 0.399 ms/frame.
4. `grd_cull_flags` — 0.1% of the highest-render window, 0.337 ms/frame.

Native quality/state evidence in the historical campaign predates the added shading features. It supplies no directly comparable native/WASM timing for this measured candidate.

## Caveats

- CYCCNT free-runs, so foreground scopes include ISR time; wake includes nested pack/submit work.
- `filter_blend` can inherit the first active parent and disappear with an inactive parent. This capture contains no such counter, and no per-pixel profiling scopes are enabled.
- Shipping selective-O3 covers GS physics, shading and hot-flag kernels plus shared lattice/orientation and renderer/driver paths. Global-O3 is a single-effect compiler comparison, not the full-roster shipping configuration.
- The extended epoch prevents reinitialization during capture. No dwell compression, simulation-speed override, or lifecycle completion is claimed.
- Both current images were captured on COM4 at source `7baf3cc430753e1a3693134150b282f212aa4e9e` with an empty source diff. The candidate includes the noise modifier speed and pigment scratch-lifetime corrections. The prior early-morning captures used COM3 and earlier artwork; their timing delta is not a controlled optimization comparison.
- Setup frame 1 is excluded from runtime rows; the entire first window is excluded from scope/ISR summaries. Complete individual frame rows in the unfinished final window remain in runtime statistics. Integer-microsecond ISR totals introduce small quantization error.
- Portable evidence manifest (supporting artifact removed) retains capture, compiler/build/environment records and original hashes. Host process environment dictionaries are removed; footer environment hashes identify the sanitized retained dumps, while measurement rows are unchanged. ELF/map artifacts remain in the provenance-named local archive.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`. the explicit `profile_o3` Reproduce command is the basic shortcut; use the Setup command for this duration and epoch under the shared-device lock.

## Global -O3 vs selective -O3

Shipping mean render is 262.576 ms/frame versus 212.279 ms/frame globally optimized (1.237×). Shipping peak is 313.695 ms versus 279.686 ms (1.122×). Both captures spill every live frame; global-O3 does not restore the 16 fps deadline. Their wall-duration runs cover different simulation-frame ranges, so these whole-capture ratios are descriptive rather than matched-frame speedups.

Global-O3 changes single-effect FLASH code by +16,024 B and ITCM by +10,832 B relative to the shipping image.

Source reachability: capture `7baf3cc4307` maps to landed `daf63d812` on a
different base (319 files differ); the captured tree is not available from the
published branch. Supporting artifacts are no longer retained.
